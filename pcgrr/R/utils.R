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
    tier_df, c("AF_TUMOR","TIER"),
    only_colnames = F, quiet = T)
  tier_df <- tier_df |>
    dplyr::select(AF_TUMOR, TIER) |>
    dplyr::mutate(TIER = paste0("TIER ", TIER))
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
    # TIER <- "NONCODING"
    # df <- data.frame(bin_name = bin_name,
    #                  bin_start = bin_start,
    #                  bin_end = bin_end,
    #                  bin = as.integer(i),
    #                  TIER = TIER, stringsAsFactors = F)
    # af_bin_df <- rbind(af_bin_df, df)
    bin_start <- bin_end
    i <- i + 1
  }

  tier_df_trans <- tier_df |>
    dplyr::mutate(
      bin = cut(.data$AF_TUMOR,
                breaks = seq(0,1,bin_size),
                right = F, include.lowest = T, labels = F))

  tier_df_trans_bin <- as.data.frame(
    dplyr::group_by(tier_df_trans, .data$TIER, .data$bin) |>
      dplyr::summarise(Count = dplyr::n(),
                       .groups = "drop"))
  af_bin_df <- af_bin_df |>
    dplyr::left_join(tier_df_trans_bin, by = c("bin", "TIER")) |>
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



#' Function that generate snv/indel + coding/noncoding + tier stats for
#' a given variant set
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
    for (n in c("n", "n_snv", "n_indel",
                "n_coding", "n_noncoding")) {
      call_stats[[name]][[n]] <- 0
    }

    if("VARIANT_CLASS" %in% colnames(callset$variant)){

      call_stats[[name]][["n"]] <- callset$variant |>
        nrow()
      call_stats[[name]][["n_snv"]] <- callset$variant |>
        dplyr::filter(.data$VARIANT_CLASS == "SNV") |>
        nrow()
      call_stats[[name]][["n_indel"]] <- callset$variant |>
        dplyr::filter(.data$VARIANT_CLASS != "SNV") |>
        nrow()
      call_stats[[name]][["n_coding"]] <- callset$variant |>
        dplyr::filter(.data$CODING_STATUS == "coding") |>
        nrow()
      call_stats[[name]][["n_noncoding"]] <- callset$variant |>
        dplyr::filter(.data$CODING_STATUS != "coding") |>
        nrow()
    }
  }

  if("TIER" %in% colnames(callset$variant)){

    for (n in c("n_tier1", "n_tier2", "n_tier3",
                "n_tier4", "n_tier5")) {
      call_stats[[name]][[n]] <- 0
    }

    call_stats[[name]][["n_tier1"]] <- callset$variant |>
      dplyr::filter(!is.na(TIER) & TIER == 1) |> nrow()
    call_stats[[name]][["n_tier2"]] <- callset$variant |>
      dplyr::filter(!is.na(TIER) & TIER == 2) |> nrow()
    call_stats[[name]][["n_tier3"]] <- callset$variant |>
      dplyr::filter(!is.na(TIER) & TIER == 3) |> nrow()
    call_stats[[name]][["n_tier4"]] <- callset$variant |>
      dplyr::filter(!is.na(TIER) & TIER == 4) |> nrow()
    call_stats[[name]][["n_tier5"]] <- callset$variant |>
      dplyr::filter(!is.na(TIER) & TIER == 5) |> nrow()
  }

  if(vartype == 'cna' &
     "VARIANT_CLASS" %in% colnames(callset$variant)){
    for (n in c("n_tsg_loss",
                "n_oncogene_gain",
                "n_other_drugtarget_gain")) {
      call_stats[[name]][[n]] <- 0
    }
    call_stats[[name]][["n_tsg_loss"]] <-
      callset$variant |>
      dplyr::filter(
        !is.na(TIER) & TIER == 3 &
          VARIANT_CLASS == "homdel") |> nrow()
    call_stats[[name]][["n_oncogene_gain"]] <-
      callset$variant |>
      dplyr::filter(
        !is.na(TIER) & TIER == 3 &
          VARIANT_CLASS == "gain") |> nrow()

    call_stats[[name]][["n_other_drugtarget_gain"]] <-
      callset$variant |>
      dplyr::filter(
        is.na(TIER) &
          VARIANT_CLASS == "gain" &
          !is.na(TARGETED_INHIBITORS_ALL2)) |>
      nrow()

  }

  return(call_stats)
}

#' Function that appends multiple HTML annotation links to variant identifiers
#' e.g. COSMIC, CLINVAR, REFSEQ etc
#'
#' @param vcf_data_df data frame with variant entries
#' @param skip elements to be ignored during annotation
#'
#' @export
append_annotation_links <- function(vcf_data_df,
                                    skip = NULL) {
  i <- 1
  while (i <= nrow(pcgrr::variant_db_url)) {

    name <- pcgrr::variant_db_url[i, ]$name
    if (!name %in% skip) {
      #pcgrr::log4r_info(paste0("Adding annotation links - ", name))
      group_by_var <- pcgrr::variant_db_url[i, ]$group_by_var
      url_prefix <- pcgrr::variant_db_url[i, ]$url_prefix
      link_key_var <- pcgrr::variant_db_url[i, ]$link_key_var
      link_display_var <- pcgrr::variant_db_url[i, ]$link_display_var
      if (!(name %in% colnames(vcf_data_df))) {
        annotation_links <-
          pcgrr::generate_annotation_link(
            vcf_data_df,
            vardb = name,
            group_by_var = group_by_var,
            url_prefix = url_prefix,
            link_key_var = link_key_var,
            link_display_var = link_display_var
          )
        if (nrow(annotation_links) > 0) {
          vcf_data_df <- vcf_data_df |>
            dplyr::left_join(dplyr::rename(annotation_links,
                                           !!rlang::sym(name) := .data$link),
                             by = c("VAR_ID"))
        }else{
          vcf_data_df[, name] <- NA
        }
      }
    }
    i <- i + 1
  }
  return(vcf_data_df)
}

#' Function that adds TCGA annotations (cohort, frequency etc.) to variant identifiers
#'
#' @param var_df data frame with variants
#' @param linktype type of link
#'
#' @return var_df
#'
#' @export
append_tcga_var_link <- function(var_df,
                                 linktype = "dbsource") {


  if (any(grepl(paste0("^TCGA_FREQUENCY$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df)))) {

    #pcgrr::log4r_info("Adding annotation links - TCGA")

    var_df_unique_slim <- dplyr::select(
      var_df, .data$VAR_ID, .data$TCGA_FREQUENCY) |>
      dplyr::filter(!is.na(.data$TCGA_FREQUENCY)) |>
      dplyr::distinct()
    if (NROW(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- var_df_unique_slim |>
        tidyr::separate_rows(.data$TCGA_FREQUENCY, sep = ",") |>
        tidyr::separate(
          .data$TCGA_FREQUENCY,
          c("tcga_cancer_code", "percentage",
            "affected", "cohort_size"),
          sep = "\\|", convert = T) |>
        dplyr::left_join(
          pcgrr::tcga_cohorts, by = "tcga_cancer_code") |>
        dplyr::arrange(
          .data$VAR_ID,
          dplyr::desc(.data$percentage))
      if (linktype == "dbsource") {
        var_df_unique_slim_melted <- var_df_unique_slim_melted |>
          dplyr::mutate(tmp_assoc = paste0(
            "<a href='https://portal.gdc.cancer.gov/projects/TCGA-",
            .data$tcga_cancer_code, "' target=\"_blank\">",
            .data$tcga_cancer_name, "</a>: ",
            .data$percentage, "% (",
            .data$affected, "/", .data$cohort_size, ")"))
      }

      var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) |>
        dplyr::summarise(TCGALINK = unlist(paste(.data$tmp_assoc, collapse = ", ")),
                         .groups = "drop") |>
        dplyr::select(.data$VAR_ID, .data$TCGALINK) |>
        magrittr::set_colnames(c("VAR_ID", "TCGA_FREQUENCY"))
      var_df <- dplyr::rename(var_df, TCGA_FREQUENCY_RAW = .data$TCGA_FREQUENCY)
      var_df <- dplyr::left_join(var_df, var_df_links,
                                 by = c("VAR_ID" = "VAR_ID"))
    }else{
      var_df$TCGA_FREQUENCY_RAW <- var_df$TCGA_FREQUENCY
    }
  }
  return(var_df)
}

#' Function that adds TFBS annotations (dbMTS) to genetic variant identifiers
#'
#' @param var_df data frame with variants
#'
#' @export
append_tfbs_annotation <-
  function(var_df) {

    if (any(grepl(paste0("^CONSEQUENCE$"), names(var_df))) &
        any(grepl(paste0("^VAR_ID$"), names(var_df))) &
        any(grepl(paste0("^REGULATORY_ANNOTATION$"), names(var_df)))) {

      #pcgrr::log4r_info("Adding TF binding site annotations - VEP regulatory")

      var_df_unique_slim <-
        dplyr::select(var_df, .data$VAR_ID,
                      .data$REGULATORY_ANNOTATION,
                      .data$CONSEQUENCE) |>
        dplyr::filter(!is.na(.data$REGULATORY_ANNOTATION) &
                        stringr::str_detect(
                          .data$CONSEQUENCE,
                          "5_prime|upstream"
                        )) |>
        dplyr::distinct()

      if (nrow(var_df_unique_slim) > 0) {
        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim |>
          tidyr::separate_rows(.data$REGULATORY_ANNOTATION, sep=",") |>
          dplyr::filter(
            stringr::str_detect(
              .data$REGULATORY_ANNOTATION, "TF_binding_site_variant"
            ))
        )

        if (nrow(var_df_unique_slim_melted) > 0) {

          pcgrr::log4r_info(paste0(
            "Found TF binding site annotations for ",
            nrow(var_df_unique_slim)," variants"))

          var_df_unique_slim_melted <- as.data.frame(
            var_df_unique_slim_melted |>
              dplyr::mutate(
                REGULATORY_ANNOTATION = stringr::str_replace(
                  .data$REGULATORY_ANNOTATION,
                  "TF_binding_site_variant\\|MotifFeature\\|ENSM0[0-9]{1,}\\|",
                  "")
              ) |>
              tidyr::separate(.data$REGULATORY_ANNOTATION,
                              into = c('cons','matrix','motif_pos',
                                       'high_inf_pos','motif_score_change',
                                       'transcription_factors'),
                              sep = "\\|",
                              remove = T) |>

              tidyr::separate_rows(.data$transcription_factors) |>
              dplyr::mutate(
                TF_BINDING_SITE_VARIANT = dplyr::case_when(
                  .data$high_inf_pos == "N" ~ "Overlap: non-critical motif position",
                  .data$high_inf_pos == "Y" ~ "Overlap: critical motif position",
                  TRUE ~ as.character(NA)
                )
              ) |>
              dplyr::mutate(
                TF_BINDING_SITE_VARIANT_INFO =
                  paste(.data$transcription_factors, .data$matrix,
                         .data$motif_pos, .data$motif_score_change,
                         .data$high_inf_pos, sep="|")
              )
          )

          var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) |>
            dplyr::summarise(
              TF_BINDING_SITE_VARIANT = paste(unique(sort(.data$TF_BINDING_SITE_VARIANT)),
                                                  collapse = ", "),
              TF_BINDING_SITE_VARIANT_INFO = paste(unique(
                .data$TF_BINDING_SITE_VARIANT_INFO),
                collapse = ", "),
              .groups = "drop") |>
            dplyr::select(.data$VAR_ID, .data$TF_BINDING_SITE_VARIANT,
                          .data$TF_BINDING_SITE_VARIANT_INFO)

          var_df <- dplyr::left_join(var_df, var_df_links,
                                     by = c("VAR_ID" = "VAR_ID"))

        }else{
          var_df$TF_BINDING_SITE_VARIANT <- NA
          var_df$TF_BINDING_SITE_VARIANT_INFO <- NA
        }
      }else{
        var_df$TF_BINDING_SITE_VARIANT <- NA
        var_df$TF_BINDING_SITE_VARIANT_INFO <- NA
      }
    }
    return(var_df)
  }

#' Function that adds miRNA target annotations (dbMTS) to
#' genetic variant identifiers
#'
#' @param var_df data frame with variants
#'
#' @export
append_dbmts_var_link <-
  function(var_df) {


    if (any(grepl(paste0("^DBMTS$"), names(var_df))) &
        any(grepl(paste0("^VAR_ID$"), names(var_df))) &
        any(grepl(paste0("^ENSEMBL_TRANSCRIPT_ID$"), names(var_df)))) {

      #pcgrr::log4r_info("Adding miRNA target site annotations (gain/loss) - dbMTS")

      var_df_unique_slim <-
        dplyr::select(
          var_df, .data$VAR_ID, .data$CLINVAR_CLASSIFICATION,
          .data$DBMTS, .data$ENSEMBL_TRANSCRIPT_ID) |>
        dplyr::filter(!is.na(.data$DBMTS) &
                        !is.na(.data$ENSEMBL_TRANSCRIPT_ID)) |>
        dplyr::distinct()
      if (nrow(var_df_unique_slim) > 0) {

        pcgrr::log4r_info(paste0(
          "Found miRNA target site annotations for ",
          nrow(var_df_unique_slim)," variants"))

        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim |>
            tidyr::separate_rows(.data$DBMTS, sep = ",") |>
            tidyr::separate(.data$DBMTS, c("ens_trans_id", "mirna_id",
                                     "algorithms", "algorithms_call",
                                     "consensus_call"),
                            sep = "\\|", convert = T) |>
            dplyr::filter(.data$ens_trans_id == .data$ENSEMBL_TRANSCRIPT_ID)
        )
        if (nrow(var_df_unique_slim_melted) > 0) {
          var_df_unique_slim_melted <- var_df_unique_slim_melted |>
            dplyr::select(-c(.data$ENSEMBL_TRANSCRIPT_ID, .data$algorithms_call)) |>
            dplyr::mutate(miRNA_TARGET_HIT = dplyr::case_when(
              .data$consensus_call == "G" ~ "gain",
              .data$consensus_call == "L" ~ "loss",
              TRUE ~ as.character(NA)
            )) |>
            dplyr::mutate(
              algorithms = stringr::str_replace_all(
                stringr::str_replace(
                  stringr::str_replace(
                    stringr::str_replace(
                      .data$algorithms, "R","RNAHybrid"),
                    "TS","TargetScan"),
                  "M","miRanda"),
                "&"," / ")
            ) |>
            dplyr::mutate(
              miRNA_TARGET_HIT_PREDICTION =
                paste0("<a href='http://www.mirbase.org/cgi-bin/mirna_entry.pl?id",
                       "=",.data$mirna_id,"' target='_blank'>",.data$mirna_id,"</a> - ", .data$miRNA_TARGET_HIT,
                       " (",.data$algorithms,")")
            )


          var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) |>
            dplyr::summarise(
              miRNA_TARGET_HIT_PREDICTION = paste(.data$miRNA_TARGET_HIT_PREDICTION,
                                                  collapse = ", "),
              miRNA_TARGET_HIT = paste(unique(.data$miRNA_TARGET_HIT),
                                       collapse = ", "),
              .groups = "drop") |>
            dplyr::select(.data$VAR_ID, .data$miRNA_TARGET_HIT, .data$miRNA_TARGET_HIT_PREDICTION)

          var_df <- dplyr::left_join(var_df, var_df_links,
                                     by = c("VAR_ID" = "VAR_ID"))
        }else{
          var_df$miRNA_TARGET_HIT_PREDICTION <- NA
          var_df$miRNA_TARGET_HIT <- NA
        }
      }else{
        var_df$miRNA_TARGET_HIT_PREDICTION <- NA
        var_df$miRNA_TARGET_HIT <- NA
      }
    }

    return(var_df)
  }

#' Function that assigns HTML links to dbNSFP prediction entries
#'
#' @param var_df data frame with variant entries
#'
#' @return var_df
#'
#' @export
append_dbnsfp_var_link <- function(var_df) {

  #pcgrr::log4r_info("Adding annotation links - dbNSFP")

  if (any(grepl(paste0("EFFECT_PREDICTIONS"), names(var_df)))) {
    var_df <- var_df |>
      dplyr::mutate(PREDICTED_EFFECT = .data$EFFECT_PREDICTIONS) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":D,", ":Damaging,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":T,", ":Tolerated,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":SN,", ":SplicingNeutral,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":AS,", ":AffectSplicing,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":PD,", ":ProbablyDamaging,"
      ))
    i <- 1
    while (i <= nrow(pcgrr::effect_prediction_algos)) {
      str_to_replace <-
        paste0(pcgrr::effect_prediction_algos[i, "algorithm"], ":")
      replacement_str <-
        paste0("<a href='",
               pcgrr::effect_prediction_algos[i, "url"], "' target='_blank'>",
               pcgrr::effect_prediction_algos[i, "display_name"], "</a>:")
      algorithm_display <-
        paste0(pcgrr::effect_prediction_algos[i, "display_name"], ":")
      var_df <- var_df |>
        dplyr::mutate(
          PREDICTED_EFFECT =
            stringr::str_replace(.data$PREDICTED_EFFECT,
                                 str_to_replace, replacement_str))
      i <- i + 1
    }
  }
  # else{
  #   var_df$PREDICTED_EFFECT <- NA
  # }
  return(var_df)

}

#' Function that adds link to targeted drugs (on and off-label)
#' for a list of variants and associated targeted
#'
#' @param vcf_data_df data frame with variants
#' @param ref_data PCGR/CPSR reference data list object
#'
#' @export
append_drug_var_link <- function(
    vcf_data_df,
    primary_site = "Breast",
    ref_data = NULL) {

  #pcgrr::log4r_info("Adding annotation links - targeted antineoplastic inhibitors")

  if (any(grepl(paste0("^SYMBOL$"), names(vcf_data_df))) &
      any(grepl(paste0("^VAR_ID$"), names(vcf_data_df))) &
      !is.null(ref_data)) {
    var_drug_df <- dplyr::select(
      vcf_data_df, c("VAR_ID","SYMBOL")) |>
      dplyr::filter(!is.na(.data$SYMBOL)) |>
      dplyr::distinct()
    if (nrow(var_drug_df) > 0) {
      var_drug_df <- var_drug_df |>
        dplyr::left_join(
          dplyr::filter(
            ref_data[['drug']][['inhibitors_any_label']],
            .data$QUERY_SITE == primary_site),
          by = c("SYMBOL")) |>
        dplyr::select(-c("QUERY_SITE")) |>
        dplyr::left_join(
          dplyr::filter(
            ref_data[['drug']][['inhibitors_on_label']],
            .data$QUERY_SITE == primary_site),
          by = c("SYMBOL")) |>
        dplyr::select(-c("QUERY_SITE")) |>
        dplyr::distinct() |>
        dplyr::filter(
          !is.na(.data$TARGETED_INHIBITORS_ALL2)) |>
        dplyr::distinct()
      if (NROW(var_drug_df) > 0) {
        vcf_data_df <-
          dplyr::left_join(
            vcf_data_df,
            var_drug_df,
            by = c("VAR_ID","SYMBOL")
          )
      }else{
        #vcf_data_df$TARGETED_INHIBITORS_OFF_LABEL <- NA
        vcf_data_df$TARGETED_INHIBITORS <- NA
        vcf_data_df$TARGETED_INHIBITORS2 <- NA
        vcf_data_df$TARGETED_INHIBITORS_ALL2 <- NA
        vcf_data_df$TARGETED_INHIBITORS_ALL <- NA
      }
    }
    else{
      vcf_data_df$TARGETED_INHIBITORS <- NA
      vcf_data_df$TARGETED_INHIBITORS2 <- NA
      vcf_data_df$TARGETED_INHIBITORS_ALL2 <- NA
      vcf_data_df$TARGETED_INHIBITORS_ALL <- NA
    }
  }

  return(vcf_data_df)
}

#' Function that appends cancer gene evidence links
#'
#' @param vcf_data_df Data frame of sample variants from VCF
#' @param ref_data PCGR reference data bundle object
#' @param primary_site Primary tumor site
#' @param pos_var variable reflecting chromosome order (POS/SEGMENT_START)
#' @return vcf_data_df
#'
#' @export
append_cancer_gene_evidence <-
  function(vcf_data_df = NULL,
           ref_data = NULL,
           primary_site = 'Any',
           pos_var = 'POS') {

    if (any(grepl(paste0("^ENTREZGENE$"), names(vcf_data_df))) &
        any(grepl(paste0("^ENSEMBL_GENE_ID$"), names(vcf_data_df)))) {

      #pcgrr::log4r_info(paste0("Adding literature evidence for cancer-relevant genes"))


      tissue_gene_ranks <- ref_data[['gene']][['otp_rank']] |>
        dplyr::select(c("ENTREZGENE", "PRIMARY_SITE", "TISSUE_ASSOC_RANK")) |>
        dplyr::filter(.data$PRIMARY_SITE == primary_site) |>
        dplyr::distinct()

      global_gene_ranks <- ref_data[['gene']][['otp_rank']] |>
        dplyr::select(c("ENTREZGENE", "GLOBAL_ASSOC_RANK")) |>
        dplyr::distinct()


      vcf_data_df_1 <- vcf_data_df |>
        dplyr::filter(!is.na(.data$ENTREZGENE))
      vcf_data_df_2 <- vcf_data_df |>
          dplyr::filter(
            is.na(.data$ENTREZGENE) &
              !is.na(.data$ENSEMBL_GENE_ID))
      vcf_data_df_3 <- vcf_data_df |>
        dplyr::filter(
          is.na(.data$ENTREZGENE) &
            is.na(.data$ENSEMBL_GENE_ID))

      if (NROW(vcf_data_df_1) > 0) {

        vcf_data_df_1 <- vcf_data_df_1 |>
          dplyr::left_join(
            dplyr::filter(
              dplyr::select(
                ref_data[["gene"]][["gene_xref"]],
                c("ENTREZGENE",
                  "ENSEMBL_GENE_ID",
                  "CANCERGENE_EVIDENCE")),
              !is.na(.data$ENTREZGENE)),
            by = c("ENTREZGENE" = "ENTREZGENE",
                   "ENSEMBL_GENE_ID" = "ENSEMBL_GENE_ID")) |>
          dplyr::distinct()

        ## Add gene ranks (Open Targets Platform)
        ## - according to primary tumor types/sites
        ## - globally (across all tumor types/sites)
        vcf_data_df_1 <- vcf_data_df_1 |>
          dplyr::left_join(
            global_gene_ranks, by = "ENTREZGENE") |>
          dplyr::mutate(GLOBAL_ASSOC_RANK = dplyr::if_else(
            is.na(.data$GLOBAL_ASSOC_RANK),
            as.numeric(0),
            as.numeric(.data$GLOBAL_ASSOC_RANK)
          ))
        if (NROW(tissue_gene_ranks) > 0) {
          tissue_gene_ranks$PRIMARY_SITE <- NULL
          vcf_data_df_1 <- vcf_data_df_1 |>
            dplyr::left_join(
              tissue_gene_ranks, by = "ENTREZGENE") |>
            dplyr::mutate(TISSUE_ASSOC_RANK = dplyr::if_else(
              is.na(.data$TISSUE_ASSOC_RANK),
              as.numeric(0),
              as.numeric(.data$TISSUE_ASSOC_RANK)
            ))
        }else{
          vcf_data_df_1 <- vcf_data_df_1 |>
            dplyr::mutate(TISSUE_ASSOC_RANK = as.numeric(0))
        }

      }

      if (NROW(vcf_data_df_2) > 0) {

        vcf_data_df_2 <- vcf_data_df_2 |>
          dplyr::left_join(
            dplyr::filter(
              dplyr::select(
                ref_data[["gene"]][["gene_xref"]],
                c("ENSEMBL_GENE_ID","CANCERGENE_EVIDENCE")),
              !is.na(.data$ENSEMBL_GENE_ID)),
            by = c("ENSEMBL_GENE_ID")) |>
          dplyr::distinct()

      }

      vcf_data_df <- dplyr::bind_rows(
        vcf_data_df_1,
        vcf_data_df_2,
        vcf_data_df_3) |>
        dplyr::distinct() |>
        dplyr::mutate(CANCERGENE_EVIDENCE = dplyr::if_else(
          .data$CANCERGENE_EVIDENCE == ".",
          as.character(NA),
          as.character(.data$CANCERGENE_EVIDENCE)
        )) |>
        pcgrr::order_variants(pos_var = pos_var)

    }

    return(vcf_data_df)

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
  if (!any(is.na(vcf_df$AF_TUMOR))) {
    vcf_df <- dplyr::filter(
      vcf_df,
      .data$AF_TUMOR >= config$somatic_snv$allelic_support$tumor_af_min)
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


#' A function that generates a HTML link for selected
#' identifiers (DBSNP, COSMIC, CLINVAR, ENTREZ)
#'
#' @param df data frame
#' @param vardb type of database
#' @param group_by_var variable used for grouping (VAR_ID)
#' @param url_prefix url prefix for link generation
#' @param link_key_var variable used in url for linking
#' @param link_display_var variable used in url for display
#' @return df_annotation_links
#'
#' @export
generate_annotation_link <- function(
  df,
  vardb = "DBSNP",
  group_by_var = "VAR_ID",
  url_prefix = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
  link_key_var = "DBSNPRSID",
  link_display_var = "DBSNPRSID") {

  df_annotation_links <- data.frame()
  invisible(assertthat::assert_that(
    is.data.frame(df), msg = paste0("Object df must be of type data.frame")
    )
  )
  assertable::assert_colnames(df, c(group_by_var, link_key_var,
                                    link_display_var),
                              only_colnames = F, quiet = T)

  if (group_by_var %in% colnames(df) &
      link_key_var %in% colnames(df) &
      link_display_var %in% colnames(df)) {
    selected_vars <- unique(c(group_by_var,
                              link_key_var,
                              link_display_var))
    tmp_df <- df[, selected_vars]
    tmp_df <- tmp_df[!is.na(tmp_df[, link_key_var]) &
                       !is.na(tmp_df[, link_display_var]), ] |>
      dplyr::distinct()


    if (nrow(tmp_df) > 0) {
      df_annotation_links <- data.frame()
      #if (vardb == "DBSNP") {
      df_annotation_links <- as.data.frame(
        tmp_df |>
          tidyr::separate_rows(!!rlang::sym(link_key_var), sep = "&|,") |>
          dplyr::mutate(
            tmp_link = paste0(
              "<a href='", url_prefix,
              !!rlang::sym(link_key_var), "' target='_blank'>",
              !!rlang::sym(link_display_var), "</a>")
          )
      )

      if (nrow(df_annotation_links) > 0) {
        df_annotation_links <- as.data.frame(
          df_annotation_links |>
            dplyr::group_by(.data$VAR_ID) |>
            dplyr::summarise(link = unlist(paste(.data$tmp_link, collapse = ", ")),
                             .groups = "drop")
        )
      }

    }
  }
  else{
    cat("WARNING: Could not generate HTML URL links", sep = "\n")
    cat("- missing url_key or grouping varable in df", sep = "\n")
  }
  return(df_annotation_links)
}


#' Function that writes a VCF intended for mutational signature analysis
#'
#' @param calls data frame with calls
#' @param sample_name sample name
#' @param output_directory Output directory for output file
#' @param vcf_fname filename for VCF
#'
#'
#' @export
write_processed_vcf <- function(calls,
                                sample_name = NULL,
                                output_directory = NULL,
                                vcf_fname = NULL) {

  pcgrr::log4r_info("Writing VCF file with input calls for signature analysis")

  header_lines <-
    c("##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE_ID,Number=1,Type=String,",
             "Description=\"Sample identifier\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

  vcf_df <- calls |>
    dplyr::select(.data$CHROM, .data$POS, .data$REF, .data$ALT) |>
    dplyr::filter(nchar(.data$REF) == 1 & nchar(.data$ALT) == 1) |>
    dplyr::filter(stringr::str_detect(.data$REF,"^(A|C|T|G)$") |
                  stringr::str_detect(.data$ALT,"^(A|C|T|G)$")) |>
    dplyr::mutate(QUAL = ".", FILTER = "PASS", ID = ".",
                  INFO = paste0("SAMPLE_ID=", sample_name)) |>
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
