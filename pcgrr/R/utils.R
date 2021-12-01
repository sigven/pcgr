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
tier_af_distribution <- function(tier_df, bin_size = 0.1) {

  af_bin_df <- data.frame()
  assertable::assert_colnames(tier_df, c("AF_TUMOR","TIER"),
                              only_colnames = F, quiet = T)
  i <- 1
  num_bins <- as.integer(1 / bin_size)
  bin_start <- 0
  while (i <= num_bins) {
    bin_end <- bin_start + bin_size
    bin_name <- as.character(paste0(bin_start, " - ", bin_end))
    j <- 1
    while (j <= 4) {
      TIER <- paste0("TIER ", j)
      df <- data.frame(bin_name = bin_name,
                       bin_start = bin_start,
                       bin_end = bin_end,
                       bin = as.integer(i),
                       TIER = TIER, stringsAsFactors = F)
      af_bin_df <- rbind(af_bin_df, df)
      j <- j + 1
    }
    TIER <- "NONCODING"
    df <- data.frame(bin_name = bin_name,
                     bin_start = bin_start,
                     bin_end = bin_end,
                     bin = as.integer(i),
                     TIER = TIER, stringsAsFactors = F)
    af_bin_df <- rbind(af_bin_df, df)
    bin_start <- bin_end
    i <- i + 1
  }

  tier_df_trans <- tier_df %>%
    dplyr::mutate(
      bin = cut(.data$AF_TUMOR,
                breaks = c(from = 0, to = 1, by = bin_size),
                right = F, include.lowest = T, labels = F))

  tier_df_trans_bin <- as.data.frame(
    dplyr::group_by(tier_df_trans, .data$TIER, .data$bin) %>%
      dplyr::summarise(Count = dplyr::n(),
                       .groups = "drop"))
  af_bin_df <- af_bin_df %>%
    dplyr::left_join(tier_df_trans_bin, by = c("bin", "TIER")) %>%
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
  vcf_df <- vcf_df %>%
    dplyr::mutate(
      !!rlang::sym(chrom_var) := as.character(!!rlang::sym(chrom_var)))
  n_before_exclusion <- nrow(vcf_df)
  nuc_chromosomes_df <- data.frame(c(as.character(seq(1:22)), "X", "Y"),
                                   stringsAsFactors = F)
  colnames(nuc_chromosomes_df) <- c(chrom_var)
  vcf_df <- dplyr::semi_join(vcf_df, nuc_chromosomes_df, by = chrom_var)
  n_after_exclusion <- nrow(vcf_df)
  log4r_info(
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
order_variants <- function(vcf_df, chrom_var = "CHROM", pos_var = "POS") {
  stopifnot(is.data.frame(vcf_df) &
              chrom_var %in% colnames(vcf_df) &
              pos_var %in% colnames(vcf_df))
  if (nrow(vcf_df) == 0)return(vcf_df)
  vcf_df %>%
    dplyr::mutate(!!rlang::sym(chrom_var) :=
                    factor(!!rlang::sym(chrom_var),
                           ordered = T,
                           levels = c(as.character(seq(1:22)), "X", "Y"))) %>%
    dplyr::arrange(!!rlang::sym(chrom_var), !!rlang::sym(pos_var)) %>%
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



#' Function that generate snv/indel + coding/noncoding stats for
#' a given variant set
#'
#' @param calls data frame with variants in predisposition_genes
#' @param name type of variant group
#'
#' @export
variant_stats_report <- function(calls, name = "v_stat") {

  call_stats <- list()
  call_stats[[name]] <- list()
  for (n in c("n", "n_snv", "n_indel", "n_coding", "n_noncoding")) {
    call_stats[[name]][[n]] <- 0
  }

  call_stats[[name]][["n"]] <- calls %>%
    nrow()
  call_stats[[name]][["n_snv"]] <- calls %>%
    dplyr::filter(.data$VARIANT_CLASS == "SNV") %>%
    nrow()
  call_stats[[name]][["n_indel"]] <- calls %>%
    dplyr::filter(.data$VARIANT_CLASS != "SNV") %>%
    nrow()
  call_stats[[name]][["n_coding"]] <- calls %>%
    dplyr::filter(.data$CODING_STATUS == "coding") %>%
    nrow()
  call_stats[[name]][["n_noncoding"]] <- calls %>%
    dplyr::filter(.data$CODING_STATUS != "coding") %>%
    nrow()

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
      log4r_info(paste0("Adding annotation links - ", name))
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
          vcf_data_df <- vcf_data_df %>%
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



#' Function that adds read support (depth, allelic fraction) for
#' tumor and normal and filters according to settings
#'
#' @param vcf_df data frame with variants
#' @param config list with workflow configuration values
#' @param precision number of significant digits for allelic fraction estimation
#'
#' @return vcf_df
#'
#' @export
append_read_support <- function(vcf_df, config = NULL, precision = 3) {

  invisible(assertthat::assert_that(
    !is.null(vcf_df), msg = "Argument 'vcf_df' cannot not be NULL"))
  invisible(assertthat::assert_that(
    is.data.frame(vcf_df),
    msg = paste0("Argument 'vcf_df' must be of type 'data.frame'")))
  invisible(assertthat::assert_that(
    !is.null(config), msg = "Argument 'config' cannot not be NULL"))
  invisible(assertthat::assert_that(
    methods::is(config, "list"),
    msg = paste0("Argument 'config' must be of type list, not ",
                 class(config)[2])))
  if (is.null(config$allelic_support))return(vcf_df)

  for (v in c("DP_TUMOR", "AF_TUMOR", "DP_CONTROL",
              "AF_CONTROL", "CALL_CONFIDENCE")) {
    vcf_df[v] <- NA
  }
  for (tag_name in names(config$allelic_support)) {
    if (config$allelic_support[[tag_name]] != "" &
        tag_name != "tumor_dp_min" &
        tag_name != "tumor_af_min" &
        tag_name != "control_dp_min" &
        tag_name != "control_af_max") {
      config$allelic_support[[tag_name]] <-
        stringr::str_replace_all(config$allelic_support[[tag_name]],
                                 "-", ".")
      if (config$allelic_support[[tag_name]] %in% colnames(vcf_df)) {
        if (tag_name == "control_af_tag") {
          vcf_df[, "AF_CONTROL"] <-
            round(as.numeric(vcf_df[, config$allelic_support[[tag_name]]]),
                  digits = precision)
        }
        if (tag_name == "control_dp_tag") {
          vcf_df[, "DP_CONTROL"] <-
            as.integer(vcf_df[, config$allelic_support[[tag_name]]])
        }
        if (tag_name == "tumor_af_tag") {
          vcf_df[, "AF_TUMOR"] <-
            round(as.numeric(vcf_df[, config$allelic_support[[tag_name]]]),
                  digits = precision)
        }
        if (tag_name == "tumor_dp_tag") {
          vcf_df[, "DP_TUMOR"] <-
            as.integer(vcf_df[, config$allelic_support[[tag_name]]])
        }
        if (tag_name == "call_conf_tag") {
          vcf_df[, "CALL_CONFIDENCE"] <-
            as.character(vcf_df[, config$allelic_support[[tag_name]]])
        }
      }
    }
  }
  return(vcf_df)
}
#' Function that appends a link to UCSC for a genomic segment
#'
#' @param var_df data frame with genomic variants
#' @param hgname name of genoome assembly ('hg38','hg19')
#' @param chrom chromosome name
#' @param start chromosome start coordinate
#' @param end chromosome end coordinate
#'
#' @return var_df
#' @export
append_ucsc_segment_link <- function(var_df,
                                     hgname = "hg38",
                                     chrom = NULL,
                                     start = NULL, end = NULL) {
  ucsc_browser_prefix <-
    paste0("http://genome.ucsc.edu/cgi-bin/hgTracks?db=", hgname, "&position=")
  if (!is.null(chrom) & !is.null(start) & !is.null(end) &
      chrom %in% colnames(var_df) & start %in% colnames(var_df) &
      end %in% colnames(var_df)) {
    var_df <- var_df %>%
      dplyr::mutate(SEGMENT_LINK = paste0(
        "<a href='", paste0(ucsc_browser_prefix,
                            paste0(!!rlang::sym(chrom),
                                   ":", !!rlang::sym(start),
                                   "-", !!rlang::sym(end)),
                            "' target=\"_blank\">",
                            paste0(!!rlang::sym(chrom), ":",
                                   !!rlang::sym(start), "-",
                                   !!rlang::sym(end)), "</a>")))
  }else{
    var_df$SEGMENT_LINK <- NA
  }
  return(var_df)


}

#' Function that adds TCGA annotations (cohort, frequency etc.) to variant identifiers
#'
#' @param var_df data frame with variants
#' @param pcgr_data PCGR data structure
#' @param linktype type of link
#'
#' @return var_df
#'
#' @export
append_tcga_var_link <- function(var_df,
                                 pcgr_data = NULL,
                                 linktype = "dbsource") {

  log4r_info("Adding annotation links - TCGA")

  if (any(grepl(paste0("^TCGA_FREQUENCY$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df))) &
      !is.null(pcgr_data)) {
    var_df_unique_slim <- dplyr::select(var_df, .data$VAR_ID, .data$TCGA_FREQUENCY) %>%
      dplyr::filter(!is.na(.data$TCGA_FREQUENCY)) %>%
      dplyr::distinct()
    if (nrow(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- var_df_unique_slim %>%
        tidyr::separate_rows(.data$TCGA_FREQUENCY, sep = ",") %>%
        tidyr::separate(.data$TCGA_FREQUENCY, c("tumor", "percentage",
                                          "affected", "cohort_size"),
                        sep = "\\|", convert = T) %>%
        dplyr::left_join(pcgr_data[["tcga"]][["projects"]], by = "tumor") %>%
        dplyr::arrange(.data$VAR_ID, dplyr::desc(.data$percentage))
      if (linktype == "dbsource") {
        var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
          dplyr::mutate(tmp_assoc = paste0(
            "<a href='https://portal.gdc.cancer.gov/projects/TCGA-",
            .data$tumor, "' target=\"_blank\">", .data$name, "</a>: ",
            .data$percentage, "% (", .data$affected, "/", .data$cohort_size, ")"))
      }

      var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) %>%
        dplyr::summarise(TCGALINK = unlist(paste(.data$tmp_assoc, collapse = ", ")),
                         .groups = "drop") %>%
        dplyr::select(.data$VAR_ID, .data$TCGALINK) %>%
        magrittr::set_colnames(c("VAR_ID", "TCGA_FREQUENCY"))
      var_df <- dplyr::rename(var_df, TCGA_FREQUENCY_RAW = .data$TCGA_FREQUENCY)
      var_df <- dplyr::left_join(var_df, var_df_links,
                                 by = c("VAR_ID" = "VAR_ID"))
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
  function(var_df){

    if (any(grepl(paste0("^CONSEQUENCE$"), names(var_df))) &
        any(grepl(paste0("^VAR_ID$"), names(var_df))) &
        any(grepl(paste0("^REGULATORY_ANNOTATION$"), names(var_df)))) {


      log4r_info("Adding TF binding site annotations for upstream and 5'UTR variants - VEP regulatory")

      var_df_unique_slim <-
        dplyr::select(var_df, .data$VAR_ID,
                      .data$REGULATORY_ANNOTATION,
                      .data$CONSEQUENCE) %>%
        dplyr::filter(!is.na(.data$REGULATORY_ANNOTATION) &
                        stringr::str_detect(
                          .data$CONSEQUENCE,
                          "5_prime|upstream"
                        )) %>%
        dplyr::distinct()

      if(nrow(var_df_unique_slim) > 0){
        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim %>%
          tidyr::separate_rows(.data$REGULATORY_ANNOTATION, sep=",") %>%
          dplyr::filter(
            stringr::str_detect(
              .data$REGULATORY_ANNOTATION, "TF_binding_site_variant"
            ))
        )

        if(nrow(var_df_unique_slim_melted) > 0){

          log4r_info(paste0(
            "Found TF binding site annotations for ",
            nrow(var_df_unique_slim)," variants"))

          var_df_unique_slim_melted <- as.data.frame(
            var_df_unique_slim_melted %>%
              dplyr::mutate(
                REGULATORY_ANNOTATION = stringr::str_replace(
                  .data$REGULATORY_ANNOTATION,
                  "TF_binding_site_variant\\|MotifFeature\\|ENSM0[0-9]{1,}\\|",
                  "")
              ) %>%
              tidyr::separate(.data$REGULATORY_ANNOTATION,
                              into = c('cons','matrix','motif_pos',
                                       'high_inf_pos','motif_score_change',
                                       'transcription_factors'),
                              sep = "\\|",
                              remove = T) %>%

              tidyr::separate_rows(.data$transcription_factors) %>%
              dplyr::mutate(
                TF_BINDING_SITE_VARIANT = dplyr::case_when(
                  .data$high_inf_pos == "N" ~ "Overlap: non-critical motif position",
                  .data$high_inf_pos == "Y" ~ "Overlap: critical motif position",
                  TRUE ~ as.character(NA)
                )
              ) %>%
              dplyr::mutate(
                TF_BINDING_SITE_VARIANT_INFO =
                  paste(.data$transcription_factors, .data$matrix,
                         .data$motif_pos, .data$motif_score_change,
                         .data$high_inf_pos, sep="|")
              )
          )

          var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) %>%
            dplyr::summarise(
              TF_BINDING_SITE_VARIANT = paste(unique(sort(.data$TF_BINDING_SITE_VARIANT)),
                                                  collapse = ", "),
              TF_BINDING_SITE_VARIANT_INFO = paste(unique(
                .data$TF_BINDING_SITE_VARIANT_INFO),
                collapse = ", "),
              .groups = "drop") %>%
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

#' Function that adds miRNA target annotations (dbMTS) to genetic variant identifiers
#'
#' @param var_df data frame with variants
#'
#' @export
append_dbmts_var_link <-
  function(var_df) {


    if (any(grepl(paste0("^DBMTS$"), names(var_df))) &
        any(grepl(paste0("^VAR_ID$"), names(var_df))) &
        any(grepl(paste0("^ENSEMBL_TRANSCRIPT_ID$"), names(var_df)))) {

      log4r_info("Adding miRNA target site annotations (gain/loss) - dbMTS")

      var_df_unique_slim <-
        dplyr::select(var_df, .data$VAR_ID, .data$CLINVAR_CLASSIFICATION,
                      .data$DBMTS, .data$ENSEMBL_TRANSCRIPT_ID) %>%
        dplyr::filter(!is.na(.data$DBMTS) & !is.na(.data$ENSEMBL_TRANSCRIPT_ID)) %>%
        dplyr::distinct()
      if (nrow(var_df_unique_slim) > 0) {

        log4r_info(paste0(
          "Found miRNA target site annotations for ",
          nrow(var_df_unique_slim)," variants"))

        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim %>%
            tidyr::separate_rows(.data$DBMTS, sep = ",") %>%
            tidyr::separate(.data$DBMTS, c("ens_trans_id", "mirna_id",
                                     "algorithms", "algorithms_call",
                                     "consensus_call"),
                            sep = "\\|", convert = T) %>%
            dplyr::filter(.data$ens_trans_id == .data$ENSEMBL_TRANSCRIPT_ID)
        )
        if(nrow(var_df_unique_slim_melted) > 0){
          var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
            dplyr::select(-c(.data$ENSEMBL_TRANSCRIPT_ID, .data$algorithms_call)) %>%
            dplyr::mutate(miRNA_TARGET_HIT = dplyr::case_when(
              .data$consensus_call == "G" ~ "gain",
              .data$consensus_call == "L" ~ "loss",
              TRUE ~ as.character(NA)
            )) %>%
            dplyr::mutate(
              algorithms = stringr::str_replace_all(
                stringr::str_replace(
                  stringr::str_replace(
                    stringr::str_replace(
                      .data$algorithms, "R","RNAHybrid"),
                    "TS","TargetScan"),
                  "M","miRanda"),
                "&"," / ")
            ) %>%
            dplyr::mutate(
              miRNA_TARGET_HIT_PREDICTION =
                paste0("<a href='http://www.mirbase.org/cgi-bin/mirna_entry.pl?id",
                       "=",.data$mirna_id,"' target='_blank'>",.data$mirna_id,"</a> - ", .data$miRNA_TARGET_HIT,
                       " (",.data$algorithms,")")
            )


          var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) %>%
            dplyr::summarise(
              miRNA_TARGET_HIT_PREDICTION = paste(.data$miRNA_TARGET_HIT_PREDICTION,
                                                  collapse = ", "),
              miRNA_TARGET_HIT = paste(unique(.data$miRNA_TARGET_HIT),
                                       collapse = ", "),
              .groups = "drop") %>%
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

  log4r_info("Adding annotation links - dbNSFP")

  if (any(grepl(paste0("EFFECT_PREDICTIONS"), names(var_df)))) {
    var_df <- var_df %>%
      dplyr::mutate(PREDICTED_EFFECT = .data$EFFECT_PREDICTIONS)
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
      var_df <- var_df %>%
        dplyr::mutate(
          PREDICTED_EFFECT =
            stringr::str_replace(.data$PREDICTED_EFFECT,
                                 str_to_replace, replacement_str))
      i <- i + 1
    }
  }
  else{
    var_df$PREDICTED_EFFECT <- NA
  }
  return(var_df)

}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param linktype type of link
#' @param pcgr_data PCGR data structure
#'
#' @export
append_drug_var_link <- function(var_df, pcgr_data = NULL,
                             linktype = "dbsource") {

  log4r_info("Adding annotation links - targeted cancer drugs")

  if (any(grepl(paste0("^CHEMBL_COMPOUND_ID$"), names(var_df))) &
      any(grepl(paste0("^SYMBOL$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df))) &
      !is.null(pcgr_data)) {
    var_df_unique_slim <- dplyr::select(var_df, .data$VAR_ID,
                                        .data$SYMBOL, .data$CHEMBL_COMPOUND_ID) %>%
      dplyr::filter(!is.na(.data$CHEMBL_COMPOUND_ID)) %>%
      dplyr::distinct()
    if (nrow(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- var_df_unique_slim %>%
        tidyr::separate_rows(.data$CHEMBL_COMPOUND_ID, sep = "&")
      chembl_drugs <-
        dplyr::select(pcgr_data[["antineopharma"]][["antineopharma"]],
                      .data$molecule_chembl_id, .data$symbol, .data$nci_concept_display_name) %>%
        dplyr::arrange(.data$symbol) %>%
        dplyr::distinct()
      var_df_unique_slim_melted <- var_df_unique_slim_melted %>%
        dplyr::left_join(chembl_drugs,
                         by = c("CHEMBL_COMPOUND_ID" = "molecule_chembl_id",
                                              "SYMBOL" = "symbol")) %>%
        dplyr::filter(!is.na(.data$nci_concept_display_name)) %>%
        dplyr::distinct()
      if (nrow(var_df_unique_slim_melted) > 0) {
        if (linktype == "dbsource") {
          var_df_unique_slim_melted <-
            var_df_unique_slim_melted %>%
            dplyr::mutate(
              tmp_antineopharma =
                paste0(
                  "<a href='https://www.targetvalidation.org/summary?drug=",
                  .data$CHEMBL_COMPOUND_ID, "' target=\"_blank\">",
                  .data$nci_concept_display_name, "</a>"))
        }
        var_df_unique_slim_melted_terms <-
          dplyr::select(var_df_unique_slim_melted,
                        .data$VAR_ID, .data$nci_concept_display_name)
        var_df_terms <- dplyr::group_by(var_df_unique_slim_melted_terms,
                                        .data$VAR_ID) %>%
          dplyr::summarise(CHEMBL_COMPOUND_TERMS =
                             paste(.data$nci_concept_display_name, collapse = ", "))
        var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) %>%
          dplyr::summarise(ANTINEOPHARMALINK =
                             unlist(paste(.data$tmp_antineopharma,
                                          collapse = ", "))) %>%
          dplyr::select(.data$VAR_ID, .data$ANTINEOPHARMALINK) %>%
          dplyr::distinct()
        var_df <- dplyr::left_join(var_df, var_df_links,
                                   by = c("VAR_ID" = "VAR_ID"))
        var_df <- dplyr::left_join(var_df, var_df_terms,
                                   by = c("VAR_ID" = "VAR_ID"))
      }else{
        var_df$ANTINEOPHARMALINK <- NA
        var_df$CHEMBL_COMPOUND_TERMS <- NA
      }
    }
    else{
      var_df$ANTINEOPHARMALINK <- NA
      var_df$CHEMBL_COMPOUND_TERMS <- NA
    }
  }
  else{
    cat(paste0("WARNING: Could not generate links with targeted compounds - ",
    "no ANTINEOPHARMA info provided in annotated VCF"), sep = "\n")
    var_df$ANTINEOPHARMALINK <- NA
    var_df$CHEMBL_COMPOUND_TERMS <- NA
  }
  return(var_df)
}

#' Function that adds HTML links to different genetic variant identifiers
#'
#' @param var_df data frame with variants
#' @param linktype type of link
#' @param pcgr_data PCGR data structure
#' @param oncotree Oncotree data frame
#'
#' @export
append_otargets_pheno_link <- function(var_df,
                                       pcgr_data = NULL,
                                       oncotree = NULL,
                                       linktype = "dbsource") {

  log4r_info(paste0("Adding annotation links - gene-cancer ",
                           "type associations (OpenTargets Platform)"))

  assertable::assert_colnames(oncotree, c("cui", "efo_id"),
                              only_colnames = F, quiet = T)

  if (any(grepl(paste0("^OPENTARGETS_DISEASE_ASSOCS$"), names(var_df))) &
      any(grepl(paste0("^ENSEMBL_GENE_ID$"), names(var_df))) &
      !is.null(pcgr_data) & nrow(oncotree) > 0) {
    var_df_unique_slim <- dplyr::select(var_df, .data$ENSEMBL_GENE_ID,
                                        .data$OPENTARGETS_DISEASE_ASSOCS) %>%
      dplyr::filter(!is.na(.data$OPENTARGETS_DISEASE_ASSOCS)) %>%
      dplyr::distinct()
    associations_found <- 0
    oncotree <- oncotree %>%
      dplyr::filter(!is.na(.data$efo_id)) %>%
      dplyr::mutate(efo_id = stringr::str_replace(.data$efo_id,":", "_"))
    if (nrow(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- as.data.frame(
        var_df_unique_slim %>%
          tidyr::separate_rows(.data$OPENTARGETS_DISEASE_ASSOCS, sep = "&") %>%
          tidyr::separate(.data$OPENTARGETS_DISEASE_ASSOCS,
                          c("efo_id", "ot_is_direct", "ot_score"),
                          sep = ":", remove = T) %>%
          dplyr::inner_join(dplyr::select(oncotree, .data$efo_id, .data$cui),
                            by = c("efo_id" = "efo_id")) %>%
          dplyr::left_join(pcgr_data[["phenotype_ontology"]][["umls"]],
                           by = c("cui" = "cui")) %>%
          dplyr::mutate(ot_score = as.numeric(.data$ot_score))
      )

      if (nrow(var_df_unique_slim_melted) > 0) {
        associations_found <- 1
        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim_melted %>%
            dplyr::group_by(.data$ENSEMBL_GENE_ID, .data$efo_id, .data$cui_name) %>%
            dplyr::summarise(score = max(.data$ot_score, na.rm = T)) %>%
            dplyr::distinct() %>%
            dplyr::arrange(dplyr::desc(.data$score))
        )

        if (linktype == "dbsource") {
          var_df_unique_slim_melted <-
            var_df_unique_slim_melted %>%
            dplyr::mutate(
              tmp_assoc =
                paste0("<a href='https://www.targetvalidation.org/evidence/",
                       .data$ENSEMBL_GENE_ID, "/", .data$efo_id, "' target=\"_blank\">",
                       .data$cui_name, "</a>"))
        }

        var_df_unique_slim_melted_terms <-
          dplyr::select(var_df_unique_slim_melted, .data$ENSEMBL_GENE_ID, .data$cui_name)
        var_df_terms <- dplyr::group_by(
          var_df_unique_slim_melted_terms, .data$ENSEMBL_GENE_ID) %>%
          dplyr::summarise(OT_DISEASE_TERMS = paste(.data$cui_name, collapse = "&"))
        var_df_links <-
          dplyr::group_by(var_df_unique_slim_melted, .data$ENSEMBL_GENE_ID) %>%
          dplyr::summarise(OT_DISEASE_LINK = unlist(paste(.data$tmp_assoc,
                                                          collapse = ", ")),
                           OPENTARGETS_RANK = max(.data$score))
        var_df_links <- dplyr::select(var_df_links,
                                      .data$ENSEMBL_GENE_ID,
                                      .data$OT_DISEASE_LINK,
                                      .data$OPENTARGETS_RANK)
        var_df <- dplyr::left_join(
          var_df, var_df_links, by = c("ENSEMBL_GENE_ID" = "ENSEMBL_GENE_ID"))
        var_df <- dplyr::left_join(
          var_df, var_df_terms, by = c("ENSEMBL_GENE_ID" = "ENSEMBL_GENE_ID"))
        var_df <- var_df %>%
          dplyr::mutate(OPENTARGETS_RANK =
                          dplyr::if_else(is.na(.data$OPENTARGETS_RANK),
                                         as.numeric(0), .data$OPENTARGETS_RANK))

      }else{
        log4r_warn(paste0("Could not generate Open Targets association links"))
        var_df$OT_DISEASE_LINK <- NA
        var_df$OT_DISEASE_TERMS <- NA
        var_df$OPENTARGETS_RANK <- 0
      }
    }else{
      if (associations_found == 0) {
        log4r_warn(paste0("Could not generate Open Targets association links"))
        var_df$OT_DISEASE_LINK <- NA
        var_df$OT_DISEASE_TERMS <- NA
        var_df$OPENTARGETS_RANK <- 0
      }
    }
  }else{
    log4r_warn(paste0("Could not generate Open Targets association ", "
                             links - no Open Targets annotations provided ", "
                             in annotated VCF"))
    var_df$OT_DISEASE_LINK <- NA
    var_df$OT_DISEASE_TERMS <- NA
    var_df$OPENTARGETS_RANK <- 0
  }
  return(var_df)
}

#' Function that adds SwissProt feature descriptions based on
#' key identifiers coming from PCGR pipeline
#'
#' @param vcf_data_df Data frame of sample variants from VCF
#' @param feature_descriptions Data frame with SwissProt feature descriptions
#' @return vcf_data_df
#'
#' @export
append_pfeature_descriptions <- function(vcf_data_df, feature_descriptions) {

  log4r_info(paste0("Extending annotation descriptions related",
                    " to UniprotKB/SwissProt protein features"))


  invisible(assertthat::assert_that(
    !is.null(vcf_data_df),
    msg = "Argument 'vcf_data_df' cannot not be NULL"))
  invisible(assertthat::assert_that(
    is.data.frame(vcf_data_df),
    msg = paste0("Argument 'vcf_data_df' must be of type 'data.frame'")))
  assertable::assert_colnames(
    vcf_data_df,
    c("UNIPROT_FEATURE", "VAR_ID", "UNIPROT_ID", "CONSEQUENCE"),
    only_colnames = F, quiet = T)
  assertable::assert_colnames(
    feature_descriptions,
    c("UNIPROT_FEATURE", "PF", "UNIPROT_ID"),
    only_colnames = T, quiet = T)

  feature_df <- dplyr::select(vcf_data_df, .data$UNIPROT_FEATURE,
                              .data$VAR_ID, .data$CONSEQUENCE, .data$UNIPROT_ID) %>%
    dplyr::filter(!is.na(.data$UNIPROT_FEATURE) & !is.na(.data$UNIPROT_ID)) %>%
    dplyr::distinct()
  if (nrow(feature_df) == 0) {
    vcf_data_df$PROTEIN_FEATURE <- NA
    return(vcf_data_df)
  }
  feature_df <- as.data.frame(
    feature_df %>%
      tidyr::separate_rows(.data$UNIPROT_FEATURE, sep = "&") %>%
      tidyr::separate_rows(.data$UNIPROT_FEATURE, sep = ",") %>%
      dplyr::left_join(
        dplyr::select(feature_descriptions, .data$UNIPROT_FEATURE, .data$UNIPROT_ID, .data$PF),
        by = c("UNIPROT_FEATURE", "UNIPROT_ID")
      ) %>%
      dplyr::filter(!is.na(.data$PF)) %>%
      dplyr::group_by(.data$VAR_ID, .data$CONSEQUENCE) %>%
      dplyr::summarise(PROTEIN_FEATURE = paste(.data$PF, collapse = ", ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        PROTEIN_FEATURE =
          dplyr::if_else(.data$PROTEIN_FEATURE == "NA",
                         as.character(NA),
                         as.character(.data$PROTEIN_FEATURE))) %>%
      dplyr::mutate(
        PROTEIN_FEATURE =
          dplyr::if_else(
            !is.na(.data$PROTEIN_FEATURE) &
              !stringr::str_detect(.data$CONSEQUENCE,
                                   "^(synonymous|missense|stop|start)"),
            as.character(NA),
            as.character(.data$PROTEIN_FEATURE))) %>%
      dplyr::filter(!is.na(.data$PROTEIN_FEATURE))
  )

  if (nrow(feature_df) > 0) {
    vcf_data_df <- dplyr::select(vcf_data_df, -.data$UNIPROT_FEATURE) %>%
      dplyr::left_join(feature_df, by = c("VAR_ID", "CONSEQUENCE"))
  }
  else{
    vcf_data_df$PROTEIN_FEATURE <- NA
  }

  return(vcf_data_df)

}

#' Function that adds GWAS citation/phenotype to GWAS hit found through PCGR annotation
#'
#' @param vcf_data_df Data frame of sample variants from VCF
#' @param gwas_citations_phenotypes data frame with variant-phenotype associations (including citations) from GWAS
#' @param p_value_threshold Required p-value to report associations from GWAS catalog
#'
#' @return vcf_data_df
#'
#' @export
append_gwas_citation_phenotype <-
  function(vcf_data_df = NULL,
           gwas_citations_phenotypes = NULL,
           p_value_threshold = 1e-6) {


    invisible(assertthat::assert_that(
      !is.null(vcf_data_df),
      msg = "Argument 'vcf_data_df' cannot not be NULL"))
    invisible(assertthat::assert_that(
      is.data.frame(vcf_data_df),
      msg = paste0("Argument 'vcf_data_df' must be of type 'data.frame'")))
    invisible(assertthat::assert_that(
      !is.null(gwas_citations_phenotypes),
      msg = "Argument 'gwas_citations_phenotypes' cannot not be NULL"))
    invisible(assertthat::assert_that(
      is.data.frame(gwas_citations_phenotypes),
      msg =
        "Argument 'gwas_citations_phenotypes' must be of type 'data.frame'"))
    assertable::assert_colnames(
      vcf_data_df, c("GWAS_HIT", "VAR_ID"),
      only_colnames = F, quiet = T)
    assertable::assert_colnames(
      gwas_citations_phenotypes,
      c("gwas_key", "GWAS_PH", "GWAS_CIT","p_value_num"),
      only_colnames = F, quiet = T)

    gwas_citations_phenotypes <- gwas_citations_phenotypes %>%
      dplyr::filter(.data$p_value_num <= .data$p_value_threshold)
    log4r_info(paste0("Adding citations/phenotypes underlying ",
                      "GWAS hits (NHGRI-EBI GWAS Catalog)"))


  if ("GWAS_HIT" %in% colnames(vcf_data_df) &
      "VAR_ID" %in% colnames(vcf_data_df)) {
    feature_df <- dplyr::select(vcf_data_df, .data$GWAS_HIT, .data$VAR_ID) %>%
      dplyr::filter(!is.na(.data$GWAS_HIT)) %>%
      dplyr::distinct()
    if (nrow(feature_df) == 0) {
      vcf_data_df$GWAS_CITATION <- NA
      vcf_data_df$GWAS_PHENOTYPE <- NA
      return(vcf_data_df)
    }


    feature_df <- as.data.frame(
      feature_df %>%
        tidyr::separate_rows(.data$GWAS_HIT, sep = ",") %>%
        tidyr::separate(.data$GWAS_HIT,
                        into = c("rsid", "risk_allele", "pmid",
                                 "tagsnp", "p_value", "efo_id"),
                        sep = "\\|", remove = F) %>%
        dplyr::mutate(gwas_key = paste(.data$rsid, .data$efo_id, .data$pmid, sep = "_")) %>%
        dplyr::left_join(dplyr::select(gwas_citations_phenotypes,
                                       .data$gwas_key, .data$GWAS_PH, .data$GWAS_CIT),
                         by = c("gwas_key")) %>%
        dplyr::filter(!is.na(.data$GWAS_CIT)) %>%
        dplyr::group_by(.data$VAR_ID) %>%
        dplyr::summarise(GWAS_PHENOTYPE = paste(unique(.data$GWAS_PH),
                                                collapse = "; "),
                         GWAS_CITATION = paste(unique(.data$GWAS_CIT),
                                               collapse = "; ")) %>%
        dplyr::mutate(GWAS_CITATION =
                        dplyr::if_else(.data$GWAS_CITATION == "NA",
                                       as.character(NA),
                                       as.character(.data$GWAS_CITATION))) %>%
        dplyr::mutate(GWAS_PHENOTYPE =
                        dplyr::if_else(.data$GWAS_PHENOTYPE == "NA",
                                       as.character(NA),
                                       as.character(.data$GWAS_PHENOTYPE))) %>%
        dplyr::filter(!is.na(.data$GWAS_CITATION) & !is.na(.data$GWAS_PHENOTYPE))
    )
    if (nrow(feature_df) > 0) {

      log4r_info(paste0(
        "Found n = ",
        nrow(feature_df),
        " variants associated with genome-wide association studies"))

      vcf_data_df <- vcf_data_df %>%
        dplyr::left_join(feature_df, by = c("VAR_ID" = "VAR_ID"))
    }else{
      vcf_data_df$GWAS_CITATION <- NA
      vcf_data_df$GWAS_PHENOTYPE <- NA
    }
  }else{
    vcf_data_df$GWAS_CITATION <- NA
    vcf_data_df$GWAS_PHENOTYPE <- NA
  }

  return(vcf_data_df)

}


#' Function that filters variant set on (depth, allelic fraction)
#' for tumor and normal and filters according to settings
#'
#' @param vcf_df data frame with variants
#' @param config list with workflow configuration values
#' @param precision number of significant digits for allelic fraction estimation
#'
#' @return vcf_df
#'
#' @export
filter_read_support <- function(vcf_df, config = NULL, precision = 3) {

  log4r_info(paste0(
    paste0("Filtering tumor variants based on allelic depth/fraction (min_dp_tumor=",
    config$allelic_support$tumor_dp_min,
    ", min_af_tumor=",
    config$allelic_support$tumor_af_min, ")"))
  )

  log4r_info(paste0(
    "Filtering tumor variants based on allelic depth/fraction (min_dp_control=",
    config$allelic_support$control_dp_min,
    ", max_af_control=",
    config$allelic_support$control_af_max, ")"))

  n_before_dp_af_filtering <- nrow(vcf_df)
  if (!any(is.na(vcf_df$DP_TUMOR))) {
    vcf_df <- dplyr::filter(vcf_df,
                            .data$DP_TUMOR >= config$allelic_support$tumor_dp_min)
  }
  if (!any(is.na(vcf_df$AF_TUMOR))) {
    vcf_df <- dplyr::filter(vcf_df,
                            .data$AF_TUMOR >= config$allelic_support$tumor_af_min)
  }
  if (!any(is.na(vcf_df$AF_CONTROL))) {
    vcf_df <- dplyr::filter(vcf_df,
                            .data$AF_CONTROL <= config$allelic_support$control_af_max)
  }
  if (!any(is.na(vcf_df$DP_CONTROL))) {
    vcf_df <- dplyr::filter(vcf_df,
                            .data$DP_CONTROL >= config$allelic_support$control_dp_min)
  }
  n_removed <- n_before_dp_af_filtering - nrow(vcf_df)
  percentage <- round(as.numeric((n_removed / n_before_dp_af_filtering) * 100),
                      digits = 2)
  log4r_info(
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
                       !is.na(tmp_df[, link_display_var]), ] %>%
      dplyr::distinct()


    if (nrow(tmp_df) > 0) {
      df_annotation_links <- data.frame()
      #if (vardb == "DBSNP") {
      df_annotation_links <- as.data.frame(
        tmp_df %>%
          tidyr::separate_rows(!!rlang::sym(link_key_var), sep = "&|,") %>%
          dplyr::mutate(
            tmp_link = paste0("<a href='", url_prefix,
                              !!rlang::sym(link_key_var), "'target=\"_blank\">",
                              !!rlang::sym(link_display_var), "</a>")
          )
      )

      if (nrow(df_annotation_links) > 0) {
        df_annotation_links <- as.data.frame(
          df_annotation_links %>%
            dplyr::group_by(.data$VAR_ID) %>%
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


#' Function that assigns explicit genotype (heterozygous/homozygous) from VCF GT tag
#'
#' @param vcf_df data frame with variant data from VCF
#'
#' @return vcf_df
#'
#' @export
determine_genotype <- function(vcf_df) {

  invisible(assertthat::assert_that(
    is.data.frame(vcf_df),
    msg = paste0("Argument 'vcf_df' must be of type data.frame, not",
                 class(vcf_df)
                 )
    )
  )
  vcf_df$GENOTYPE <- "NA"
  if ("GT" %in% colnames(vcf_df)) {
    vcf_df <- vcf_df %>%
      dplyr::mutate(GENOTYPE = dplyr::case_when(
        GT %in% pcgrr::heterozygous_states ~ "heterozygous",
        GT %in% pcgrr::homozygous_states ~ "homozygous",
        TRUE ~ as.character("ND")

      ))
  }
  return(vcf_df)
}

#' Function that performs a sanity check of query VCF for pcgrr package,
#' ensuring that all required tags are present in annotated VCF
#'
#' @param vcf_df data frame with variants from VCF
#' @param pcgr_data PCGR data object
#' @param config PCGR/CPSR config object
#' @param workflow type of workflow (pcgr/cpsr)
#'
#'
#' @export
data_integrity_check <- function(vcf_df,
                                 pcgr_data,
                                 config = NULL,
                                 workflow = "pcgr") {
  log4r_info("Verifying data integrity of input callset")

  stopifnot(is.data.frame(vcf_df) & !is.null(pcgr_data))
  stopifnot(!is.null(pcgr_data[["annotation_tags"]][["vcf_cpsr"]]) &
              !is.null(pcgr_data[["annotation_tags"]][["vcf_pcgr"]]))

  vars_required <- pcgr_data[["annotation_tags"]][["vcf_cpsr"]]
  if (workflow == "pcgr") {
    vars_required <- pcgr_data[["annotation_tags"]][["vcf_pcgr"]]
    vars_required <- vars_required[!vars_required$tag == "PANEL_OF_NORMALS", ]

    if(!is.null(config)){
      if(config$other$vep_regulatory == F){
        vars_required <- vars_required[!vars_required$tag == "REGULATORY_ANNOTATION", ]
      }
    }
  }

  vcf_required_vars <-
    dplyr::bind_rows(data.frame(tag = "CHROM", stringsAsFactors = F),
                     data.frame(tag = "REF", stringsAsFactors = F),
                     data.frame(tag = "ALT", stringsAsFactors = F),
                     data.frame(tag = "POS", stringsAsFactors = F),
                     data.frame(tag = "QUAL", stringsAsFactors = F),
                     data.frame(tag = "FILTER", stringsAsFactors = F))
  vcf_required_vars$tag <- as.character(vcf_required_vars$tag)
  vars_required <- dplyr::bind_rows(vars_required, vcf_required_vars)
  i <- 1
  while (i <= nrow(vars_required)) {
    tag <- vars_required[i, "tag"]
    if (!(tag %in% colnames(vcf_df))) {
      if (workflow == "pcgr") {
        stop(
          paste0("Missing required variable (",
                 tag, ") in annotated TSV file from PCGR workflow - quitting"))
      }else{
        stop(
          paste0("Missing required variable (",
                 tag, ") in annotated TSV file from CPSR workflow - quitting"))
      }
    }
    i <- i + 1
  }



  return(vcf_df)
}

#' Function that assigns proper logical values to logical entries
#' read from VCF file (TRUE = True (VCF), FALSE = False (VCF))
#'
#' @param vcf_data_df variant data frame
#'
#' @export
recode_logical_vars <- function(vcf_data_df) {
  for (v in c("ONCOGENE", "TUMOR_SUPPRESSOR", "NETWORK_CG",
              "LAST_EXON", "LAST_INTRON", "PANEL_OF_NORMALS",
              "NULL_VARIANT", "SPLICE_DONOR_RELEVANT",
              "WINMASKER_HIT", "SIMPLEREPEATS_HIT")) {
    if (v %in% colnames(vcf_data_df)) {
      vcf_data_df[, v] <-
        as.logical(dplyr::recode(vcf_data_df[, v],
                                 True = TRUE, False = FALSE))
    }
  }
  return(vcf_data_df)

}


#' Function that reads a fully annotated VCF from PCGR VEP/vcfanno pipeline
#'
#' @param tsv_gz_file Bgzipped VCF file
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param config Object with PCGR configuration parameters
#' @param oncotree data frame with sample-dependent phenotype terms from UMLS
#' @param cpsr logical indicating call retrieval from CPSR workflow
#' @param n_lines_skip number of lines to skip from vcf2tsv file
#' @param maf_filenames list with filenames for sample MAFs (temp and final)
#'
#' @return vcf_data_df
#'
#' @export
get_calls <- function(tsv_gz_file,
                      pcgr_data,
                      sample_name,
                      config,
                      oncotree = NULL,
                      cpsr = F,
                      n_lines_skip = 1,
                      maf_filenames = NULL) {

  ## check that arguments are valid
  invisible(assertthat::assert_that(
    !is.null(sample_name), msg = "Argument 'sample_name' cannot not be NULL"))
  invisible(assertthat::assert_that(
    methods::is(sample_name, "character"),
    msg = paste0("Argument 'sample_name' must be of type 'character', not ",
                 class(sample_name))))
  invisible(assertthat::assert_that(
    !is.null(tsv_gz_file), msg = "Argument 'tsv_gz_file' cannot not be NULL"))
  invisible(assertthat::assert_that(
    file.exists(tsv_gz_file),
    msg = paste0("File ", tsv_gz_file, " does not exist")))
  # invisible(assertthat::assert_that(
  #   nat.utils::is.gzip(tsv_gz_file) == T,
  #   msg = paste0("File ", tsv_gz_file, " is not a gzipped TSV file")))
  invisible(assertthat::assert_that(
    !is.null(pcgr_data), msg = "Argument 'pcgr_data' cannot not be NULL"))
  invisible(assertthat::assert_that(
    methods::is(pcgr_data, "list"),
    msg = paste0("Argument 'pcgr_data' must be of type list, not ",
                 class(pcgr_data))))
  invisible(
    assertthat::assert_that(!is.null(config),
                            msg = "Argument 'config' cannot not be NULL"))
  invisible(
    assertthat::assert_that(
      methods::is(config, "list"),
      msg = paste0("Argument 'config' must be of type list, not ",
                   class(config)[2])))
  invisible(assertthat::assert_that(
    !is.null(oncotree), msg = "Argument 'oncotree' cannot not be NULL"))

  ## read VCF data
  vcf_data_df <- utils::read.table(gzfile(tsv_gz_file), skip = n_lines_skip,
                            sep = "\t", header = T, stringsAsFactors = F,
                            quote = "", comment.char = "",
                            na.strings = c("."))
  if (nrow(vcf_data_df) == 0) return(vcf_data_df)

  ## Fix erroneous parsing of cases where all
  ## reference/alternate alleles are "T",
  ## converted to logical with read.table
  if (is.logical(vcf_data_df$REF)) {
    vcf_data_df$REF <- as.character("T")
  }
  if (is.logical(vcf_data_df$ALT)) {
    vcf_data_df$ALT <- as.character("T")
  }


  wflow <- "pcgr"
  if (cpsr == T) {
    wflow <- "cpsr"
  }

  ## Annotation variables added by get_calls()
  pcgr_cpsr_columns <-
    c("GENOME_VERSION",
      "PROTEIN_CHANGE",
      "CONSEQUENCE",
      "GENOMIC_CHANGE",
      "VAR_ID",
      "OPENTARGETS_ASSOCIATIONS",
      "DOCM_DISEASE",
      "DOCM_LITERATURE",
      "CLINVAR",
      "CLINVAR_TRAITS_ALL",
      "GENE_NAME",
      "GENENAME",
      "CANCERGENE_SUPPORT",
      "OPENTARGETS_RANK",
      "TARGETED_DRUGS",
      "CANCER_ASSOCIATIONS",
      "DBSNP",
      "COSMIC",
      "PROTEIN_DOMAIN",
      "CLINVAR_PHENOTYPE",
      "NCBI_REFSEQ",
      "AF_TUMOR",
      "DP_TUMOR",
      "AF_CONTROL",
      "DP_CONTROL",
      "CALL_CONFIDENCE",
      "PFAM_DOMAIN_NAME",
      "GWAS_CITATION",
      "GWAS_PHENOTYPE",
      "TF_BINDING_SITE_VARIANT",
      "TF_BINDING_SITE_VARIANT_INFO",
      "miRNA_TARGET_HIT",
      "miRNA_TARGET_HIT_PREDICTION"
    )
  vcf_data_df <- vcf_data_df[, !(colnames(vcf_data_df) %in% pcgr_cpsr_columns)]

  ## Perform data integrity check of input
  vcf_data_df <- vcf_data_df %>%
    pcgrr::data_integrity_check(pcgr_data = pcgr_data,
                                config = config,
                                workflow = wflow) %>%
    dplyr::mutate(GENOMIC_CHANGE = paste0(.data$CHROM, ":g.", .data$POS, .data$REF, ">", .data$ALT)) %>%
    dplyr::mutate(VAR_ID = paste(.data$CHROM, .data$POS, .data$REF, .data$ALT, sep = "_")) %>%
    pcgrr::detect_vcf_sample_name(cpsr = cpsr, sample_name = sample_name) %>%
    pcgrr::get_ordinary_chromosomes(chrom_var = "CHROM")
  if (nrow(vcf_data_df) == 0) return(vcf_data_df)

  vcf_data_df <- vcf_data_df %>%
    pcgrr::order_variants() %>%
    dplyr::rename(CONSEQUENCE = .data$Consequence)

  ##
  if (cpsr == T) {
    if ("LoF" %in% colnames(vcf_data_df)) {
      vcf_data_df <- vcf_data_df %>%
        dplyr::mutate(
          LOSS_OF_FUNCTION = dplyr::if_else(
            !is.na(.data$LoF) & .data$LoF == "HC",
            TRUE, FALSE, FALSE)) %>%
        ## Ignore LoF predictions for missense variants (bug in LofTee?)
        dplyr::mutate(
          LOSS_OF_FUNCTION = dplyr::if_else(
            !is.na(.data$CONSEQUENCE) &
              .data$CONSEQUENCE == "missense_variant" &
              .data$LOSS_OF_FUNCTION == T, FALSE, .data$LOSS_OF_FUNCTION))
    }
    ## Add GWAS annotations (phenotypes)
    if (!is.null(config[["gwas"]][["p_value_min"]])) {
      vcf_data_df <- vcf_data_df %>%
        pcgrr::append_gwas_citation_phenotype(
          gwas_citations_phenotypes =
            pcgr_data[["gwas"]][["citations_phenotypes"]],
          p_value_threshold = config[["gwas"]][["p_value_min"]]
          )
    }
    ## adding miRNA target site overlap (dbMTS) & TF binding site variants (VEP)
    vcf_data_df <- vcf_data_df %>%
      pcgrr::determine_genotype() %>%
      pcgrr::append_dbmts_var_link() %>%
      pcgrr::append_tfbs_annotation()

  }else{
    if (!("PANEL_OF_NORMALS" %in% colnames(vcf_data_df))) {
      vcf_data_df <- vcf_data_df %>%
        dplyr::mutate(PANEL_OF_NORMALS = "False")
    }
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_read_support(config = config)

    if ("maf_tmp" %in% names(maf_filenames) & "maf" %in% names(maf_filenames)) {
      pcgrr::update_maf_allelic_support(
        vcf_data_df, maf_filenames[["maf_tmp"]], maf_filenames[["maf"]])
    }

   vcf_data_df <- vcf_data_df %>%
      pcgrr::filter_read_support(config = config) %>%
      dplyr::mutate(
        PUTATIVE_DRIVER_MUTATION =
          dplyr::if_else(!is.na(.data$PUTATIVE_DRIVER_MUTATION),
                         TRUE, FALSE))
  }

  if (nrow(vcf_data_df) == 0) return(vcf_data_df)

  vcf_data_df <- vcf_data_df %>%
    dplyr::mutate(SYMBOL = .data$SYMBOL_ENTREZ) %>%
    dplyr::mutate(
      CLINVAR_CONFLICTED = dplyr::case_when(
        .data$CLINVAR_CONFLICTED == "1" ~ TRUE,
        .data$CLINVAR_CONFLICTED == "0" ~ FALSE,
        FALSE ~ as.logical(.data$CLINVAR_CONFLICTED))) %>%
    dplyr::mutate(
      GENOME_VERSION = pcgr_data[["assembly"]][["grch_name"]],
      PROTEIN_CHANGE = .data$HGVSp) %>%
    dplyr::mutate(
      PROTEIN_CHANGE = dplyr::if_else(
        stringr::str_detect(.data$PROTEIN_CHANGE, ":"),
        stringr::str_split_fixed(.data$PROTEIN_CHANGE, pattern = ":", 2)[, 2],
        as.character(.data$PROTEIN_CHANGE))) %>%
    dplyr::mutate(
      PROTEIN_CHANGE = dplyr::if_else(
        stringr::str_detect(.data$PROTEIN_CHANGE, "^ENSP"),
        as.character(NA),
        as.character(.data$PROTEIN_CHANGE))) %>%
    dplyr::mutate(CLINVAR_MSID = as.integer(.data$CLINVAR_MSID)) %>%
    dplyr::mutate(EXON = stringr::str_replace(.data$EXON, "-[0-9]{1,}", "")) %>%
    dplyr::mutate(
      EXON = dplyr::if_else(
        !is.na(.data$EXON) & stringr::str_detect(.data$EXON, "^([0-9]{1,}/[0-9]{1,})$"),
        as.integer(stringr::str_split_fixed(.data$EXON, "/", 2)[, 1]),
        as.integer(NA))) %>%
    pcgrr::append_pfeature_descriptions(
      feature_descriptions =
        pcgr_data[["protein_features"]][["swissprot"]]) %>%
    dplyr::mutate(
      PFAM_DOMAIN = as.character(
        stringr::str_replace(.data$PFAM_DOMAIN, "Pfam_domain:", ""))) %>%
    dplyr::left_join(
      dplyr::select(pcgr_data[["protein_domains"]][["pfam"]],
                    .data$pfam_id, .data$pfam_name),
      by = c("PFAM_DOMAIN" = "pfam_id")) %>%
    dplyr::rename(PFAM_DOMAIN_NAME = .data$pfam_name) %>%
    dplyr::left_join(pcgr_data[["biomarkers"]][["docm"]], by = c("VAR_ID")) %>%
    dplyr::mutate(ENTREZ_ID = as.character(.data$ENTREZ_ID)) %>%
    dplyr::mutate(Gene = as.character(.data$Gene)) %>%
    pcgrr::recode_logical_vars()



  log4r_info(paste0("Number of PASS variants: ", nrow(vcf_data_df)))
  if (any(grepl(paste0("VARIANT_CLASS$"), names(vcf_data_df)))) {
    n_snvs <-
      vcf_data_df %>% dplyr::filter(!is.na(.data$VARIANT_CLASS) &
                                      .data$VARIANT_CLASS == "SNV") %>%
      nrow()
    n_deletions <- vcf_data_df %>%
      dplyr::filter(
        !is.na(.data$VARIANT_CLASS) &
          (.data$VARIANT_CLASS == "deletion" | .data$VARIANT_CLASS == "indel")) %>%
      nrow()
    n_insertions <- vcf_data_df %>%
      dplyr::filter(!is.na(.data$VARIANT_CLASS) & .data$VARIANT_CLASS == "insertion") %>%
      nrow()
    n_substitutions <- vcf_data_df %>%
      dplyr::filter(!is.na(.data$VARIANT_CLASS) & .data$VARIANT_CLASS == "substitution") %>%
      nrow()
    log4r_info(
      paste0("Number of SNVs: ", n_snvs))
    log4r_info(
      paste0("Number of deletions: ", n_deletions))
    log4r_info(
      paste0("Number of insertions: ", n_insertions))
    log4r_info(
      paste0("Number of block substitutions: ", n_substitutions))
  }



  if (nrow(vcf_data_df) == 0)return(vcf_data_df)

  vcf_data_df_1 <-
    dplyr::left_join(
      dplyr::filter(vcf_data_df, !is.na(.data$ENTREZ_ID)),
      dplyr::filter(dplyr::select(pcgr_data[["gene_xref"]][["gencode"]],
                                  .data$ENTREZ_ID, .data$Gene,
                                  .data$GENENAME, .data$CANCERGENE_SUPPORT),
                                   !is.na(.data$ENTREZ_ID)),
      by = c("ENTREZ_ID", "Gene"))
  vcf_data_df_2 <-
    dplyr::left_join(
      dplyr::filter(vcf_data_df, is.na(.data$ENTREZ_ID)),
      dplyr::select(pcgr_data[["gene_xref"]][["gencode"]],
                    .data$Gene, .data$GENENAME, .data$CANCERGENE_SUPPORT),
      by = c("Gene"))
  vcf_data_df <- rbind(vcf_data_df_1, vcf_data_df_2) %>%
    pcgrr::order_variants()

  log4r_info("Extending annotation descriptions related to KEGG pathways")
  vcf_data_df <- dplyr::left_join(vcf_data_df,
                                  pcgr_data[["kegg"]][["pathway_links"]],
                                  by = c("ENTREZ_ID" = "gene_id")) %>%
    dplyr::rename(KEGG_PATHWAY = .data$kegg_pathway_urls)

  clinvar <- dplyr::select(pcgr_data[["clinvar"]][["variants"]],
                           .data$CLINVAR_TRAITS_ALL, .data$CLINVAR_MSID, .data$VAR_ID) %>%
    dplyr::filter(!is.na(.data$CLINVAR_MSID) & !is.na(.data$VAR_ID))
  if ("CLINVAR_MSID" %in% colnames(vcf_data_df) &
      "VAR_ID" %in% colnames(vcf_data_df)) {
    log4r_info("Extending annotation descriptions related to ClinVar")
    vcf_data_df <-
      dplyr::left_join(vcf_data_df, clinvar,
                       by = c("CLINVAR_MSID", "VAR_ID")) %>%
      dplyr::mutate(CLINVAR_MSID = as.character(.data$CLINVAR_MSID))
  }

  vcf_data_df <-
    pcgrr::df_string_replace(
      vcf_data_df, strings = c("CONSEQUENCE", "REFSEQ_MRNA"),
      pattern = "&", replacement = ", ",
      replace_all = T)
  vcf_data_df <-
    pcgrr::df_string_replace(
      vcf_data_df,
      strings = c("VEP_ALL_CSQ",
                  "REGULATORY_ANNOTATION",
                  "DOCM_DISEASE",
                  "MUTATION_HOTSPOT_CANCERTYPE",
                  "ICGC_PCAWG_OCCURRENCE"),
      pattern = ",", replacement = ", ", replace_all = T)


  if ("EFFECT_PREDICTIONS" %in% colnames(vcf_data_df)) {
    vcf_data_df <- vcf_data_df %>%
      dplyr::mutate(
        EFFECT_PREDICTIONS = stringr::str_replace_all(
          .data$EFFECT_PREDICTIONS, "\\.&|\\.$", "NA&"))  %>%
      dplyr::mutate(
        EFFECT_PREDICTIONS = stringr::str_replace_all(
          .data$EFFECT_PREDICTIONS, "&$", "")) %>%
      dplyr::mutate(
        EFFECT_PREDICTIONS = stringr::str_replace_all(
          .data$EFFECT_PREDICTIONS, "&", ", ")) %>%
      pcgrr::append_dbnsfp_var_link()
  }

  if ("TCGA_FREQUENCY" %in% colnames(vcf_data_df)) {
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_tcga_var_link(pcgr_data = pcgr_data)
  }

  vcf_data_df <- vcf_data_df %>%
    pcgrr::append_annotation_links()

  if (!("TARGETED_DRUGS" %in% colnames(vcf_data_df))) {
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_drug_var_link(pcgr_data = pcgr_data) %>%
      dplyr::rename(TARGETED_DRUGS = .data$ANTINEOPHARMALINK)
  }

  if (!("OPENTARGETS_ASSOCIATIONS" %in% colnames(vcf_data_df))) {
    vcf_data_df <- vcf_data_df %>%
      pcgrr::append_otargets_pheno_link(pcgr_data = pcgr_data,
                                        oncotree = oncotree) %>%
      dplyr::rename(OPENTARGETS_ASSOCIATIONS = .data$OT_DISEASE_LINK)
  }


  return(vcf_data_df)

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
write_processed_vcf <- function(calls, sample_name = NULL,
                                output_directory = NULL,
                                vcf_fname = NULL) {

  log4r_info("Writing VCF file with input calls for signature analysis")

  header_lines <-
    c("##fileformat=VCFv4.2",
      paste0("##INFO=<ID=SAMPLE_ID,Number=1,Type=String,",
             "Description=\"Sample identifier\">"),
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

  vcf_df <- calls %>%
    dplyr::select(.data$CHROM, .data$POS, .data$REF, .data$ALT) %>%
    dplyr::filter(nchar(.data$REF) == 1 & nchar(.data$ALT) == 1) %>%
    dplyr::filter(stringr::str_detect(.data$REF,"^(A|C|T|G)$") |
                  stringr::str_detect(.data$ALT,"^(A|C|T|G)$")) %>%
    dplyr::mutate(QUAL = ".", FILTER = "PASS", ID = ".",
                  INFO = paste0("SAMPLE_ID=", sample_name)) %>%
    dplyr::distinct()

  sample_vcf_content_fname <-
    file.path(output_directory,
              paste0(sample_name, ".",
                     sample(100000, 1), ".vcf_content.tsv"))
  write(header_lines, file = vcf_fname, sep = "\n")

  sample_vcf <- vcf_df[, c("CHROM", "POS", "ID", "REF",
                          "ALT", "QUAL", "FILTER", "INFO")] %>%
    dplyr::mutate(REF = as.character(.data$REF)) %>%
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
    log4r_info(paste0("Found the following VCF sample names: ",
                             paste(unique_sample_names, collapse = ", ")))

    if (length(unique_sample_names) > 1 & cpsr == T) {
      log4r_info(paste0("Found more than one sample name - VCF with somatic ",
                     "calls? Expecting single sample germline VCF for CPSR"))
      stop()
    }
  }
  df <- df %>%
    dplyr::mutate(VCF_SAMPLE_ID = sample_name)
  return(df)
}

#' Function that updates a MAF file (produced by vcf2maf) with allelic support data
#'
#' @param calls data frame with variant calls
#' @param maf_fname_tmp Filename for MAF produced by vcf2maf
#' @param maf_fname Filename for MAF that includes variant allelic support (final MAF)
#' @param delete_raw logical indicating if initial MAF should be deleted

#' @export
#'
update_maf_allelic_support <- function(calls, maf_fname_tmp,
                                       maf_fname, delete_raw = T) {

  assertable::assert_colnames(calls, c("CHROM",
                                       "POS",
                                      "DP_TUMOR",
                                      "AF_TUMOR",
                                      "DP_CONTROL",
                                      "AF_CONTROL",
                                      "VARIANT_CLASS"),
                              only_colnames = F, quiet = T)
  log4r_info(paste0("Updating MAF file with information regarding ",
                    "variant sequencing depths and allelic support"))
  maf_data <- NULL
  if (file.exists(maf_fname_tmp) & file.info(maf_fname_tmp)$size > 0) {
    maf_data <- utils::read.table(maf_fname_tmp, skip = 1, header = T,
                           sep = "\t", na.strings = c(""),
                           stringsAsFactors = F,
                           quote = "", comment.char = "") %>%
      dplyr::mutate(Chromosome = as.character(.data$Chromosome))

  allele_columns <- c("Reference_Allele",
                      "Tumor_Seq_Allele1",
                      "Tumor_Seq_Allele2",
                     "Match_Norm_Seq_Allele1",
                     "Match_Norm_Seq_Allele2",
                     "Tumor_Validation_Allele1",
                     "Tumor_Validation_Allele2",
                     "Match_Norm_Validation_Allele1",
                     "Match_Norm_Validation_Allele2")
  for (c in allele_columns) {
    if (!is.null(maf_data[, c])) {
      if (is.logical(maf_data[, c])) {
        maf_data[, c] <- as.character("T")
      }
    }
  }

  calls_maf <- calls %>%
    dplyr::select(.data$CHROM, .data$POS, .data$DP_TUMOR, .data$AF_TUMOR,
                  .data$DP_CONTROL, .data$AF_CONTROL, .data$VARIANT_CLASS) %>%
    dplyr::rename(Chromosome = .data$CHROM, Start_Position = .data$POS) %>%
    dplyr::mutate(VARIANT_CLASS = as.character(.data$VARIANT_CLASS)) %>%
    dplyr::mutate(Chromosome = as.character(.data$Chromosome)) %>%
    dplyr::mutate(Start_Position = dplyr::if_else(
      .data$VARIANT_CLASS == "deletion",
      .data$Start_Position + 1,
      as.double(.data$Start_Position)))

    if (!is.null(maf_data)) {
      if (!any(is.na(calls_maf$DP_TUMOR)) & !any(is.na(calls_maf$AF_TUMOR))) {
        calls_maf$t_depth_estimate <- calls_maf$DP_TUMOR
        calls_maf$t_ref_count_estimate <-
          calls_maf$DP_TUMOR - round(calls_maf$AF_TUMOR * calls_maf$DP_TUMOR,
                                     digits = 0)
        calls_maf$t_alt_count_estimate <-
          calls_maf$DP_TUMOR - calls_maf$t_ref_count_estimate
        calls_maf <- calls_maf %>%
          dplyr::select(-c(.data$DP_TUMOR, .data$AF_TUMOR)) %>%
          dplyr::mutate(Chromosome = as.character(.data$Chromosome))

        maf_data <- maf_data %>%
          dplyr::mutate(Chromosome = as.character(.data$Chromosome)) %>%
          dplyr::mutate(VARIANT_CLASS = as.character(.data$VARIANT_CLASS)) %>%
          dplyr::left_join(
            dplyr::select(calls_maf, -c(.data$DP_CONTROL, .data$AF_CONTROL)),
            by = c("Chromosome", "Start_Position", "VARIANT_CLASS")) %>%
          dplyr::mutate(t_depth = .data$t_depth_estimate,
                        t_ref_count = .data$t_ref_count_estimate,
                        t_alt_count = .data$t_alt_count_estimate) %>%
          dplyr::mutate(t_depth_estimate = NULL,
                        t_ref_count_estimate = NULL,
                        t_alt_count_estimate = NULL)

      }
      if (!any(is.na(calls_maf$DP_CONTROL)) &
          !any(is.na(calls_maf$AF_CONTROL))) {
        calls_maf$n_depth_estimate <- calls_maf$DP_CONTROL
        calls_maf$n_ref_count_estimate <-
          calls_maf$DP_CONTROL - round(
            calls_maf$AF_CONTROL * calls_maf$DP_CONTROL, digits = 0)
        calls_maf$n_alt_count_estimate <-
          calls_maf$DP_CONTROL - calls_maf$n_ref_count_estimate
        calls_maf <- calls_maf %>%
          dplyr::select(-c(.data$DP_CONTROL, .data$AF_CONTROL)) %>%
          dplyr::mutate(Chromosome = as.character(.data$Chromosome))


        maf_data <- maf_data %>%
          dplyr::left_join(calls_maf, by =
                             c("Chromosome",
                               "Start_Position",
                               "VARIANT_CLASS")) %>%
          dplyr::mutate(n_depth = .data$n_depth_estimate,
                        n_ref_count = .data$n_ref_count_estimate,
                        n_alt_count = .data$n_alt_count_estimate) %>%
          dplyr::mutate(n_depth_estimate = NULL,
                        n_ref_count_estimate = NULL,
                        n_alt_count_estimate = NULL)

      }

      write("#version 2.4", file = maf_fname, sep = "\n")
      utils::write.table(maf_data, file = maf_fname, col.names = T, append = T,
                  row.names = F, na = "", quote = F, sep = "\t")
      if (delete_raw == T) {
        system(paste0("rm -f ", maf_fname_tmp))
      }

    }
  }

}

#' Function that retrieves targeted drugs for a given tumor type/primary site
#'
#' @param ttype primary site/tumor type
#' @param pcgr_data PCGR data object
#' @param ignore_on_label_early_phase ignore early phase drugs (on label)
#'
#'
#' @export
targeted_drugs_pr_ttype <- function(ttype,
                                    pcgr_data,
                                    ignore_on_label_early_phase = T){

  log4r_info(
    paste0("Retrieving targeted drugs (on-label and off-label) for ",
           "indications relevant for tumor type - ", ttype))
  assertthat::assert_that(!is.null(ttype))

  if (ttype == "Cancer, NOS") {
    ttype <- "Any"
  }

  assertthat::assert_that(
    !is.null(pcgr_data[["antineopharma"]][["drug_pr_site"]][[ttype]]))

  site_candidates <- pcgr_data[["antineopharma"]][["drug_pr_site"]][[ttype]]

  drug_candidates <- list()
  drug_candidates[["on_label"]] <-
    site_candidates[["on_label"]][["late_phase"]]
  drug_candidates[["on_label_early_phase"]] <-
    site_candidates[["on_label"]][["early_phase"]]
  drug_candidates[["off_label"]] <-
    site_candidates[["off_label"]]

  if (nrow(drug_candidates[["on_label"]]) > 0) {

    drug_candidates[["on_label"]] <-
      pcgrr::targeted_drugs_summarise(
        candidate_drugs = drug_candidates[["on_label"]],
        link_label = "DRUGS_ON_LABEL",
        indication_label = "DRUGS_ON_LABEL_INDICATIONS")
  }
  if (nrow(drug_candidates[["on_label_early_phase"]]) > 0) {

    drug_candidates[["on_label_early_phase"]] <-
      pcgrr::targeted_drugs_summarise(
        candidate_drugs = drug_candidates[["on_label_early_phase"]],
        link_label = "DRUGS_ON_LABEL_EARLY_PHASE",
        indication_label = "DRUGS_ON_LABEL_EARLY_PHASE_INDICATIONS")
  }

  if (nrow(drug_candidates[["off_label"]]) > 0) {
    drug_candidates[["off_label"]] <-
      pcgrr::targeted_drugs_summarise(
        candidate_drugs = drug_candidates[["off_label"]],
        link_label = "DRUGS_OFF_LABEL",
        indication_label = "DRUGS_OFF_LABEL_INDICATIONS")
  }

  all_candidates <- dplyr::full_join(drug_candidates[["off_label"]],
                                     drug_candidates[["on_label_early_phase"]],
                                     by = "symbol")
  if (nrow(drug_candidates[["on_label"]]) > 0) {
    all_candidates <- dplyr::full_join(all_candidates,
                                       drug_candidates[["on_label"]],
                                       by = "symbol")
  }else{
    all_candidates$DRUGS_ON_LABEL <- NA
    all_candidates$DRUGS_ON_LABEL_INDICATIONS <- NA
  }

  all_drug_target_candidates <- all_candidates %>%
    dplyr::rename(SYMBOL = .data$symbol) %>%
    dplyr::arrange(.data$SYMBOL)

  if (ignore_on_label_early_phase == T) {
    all_drug_target_candidates <- all_drug_target_candidates %>%
      pcgrr::remove_cols_from_df(
        cnames = c("DRUGS_ON_LABEL_EARLY_PHASE",
                   "DRUGS_ON_LABEL_EARLY_PHASE_INDICATIONS")) %>%
      dplyr::filter(!is.na(.data$DRUGS_ON_LABEL) | !is.na(.data$DRUGS_OFF_LABEL))
  }

  return(all_drug_target_candidates)
}

#' Function that summarises available targeted drugs for a given target
#'
#' @param candidate_drugs data frame with candidate drugs
#' @param link_label column name for list of aggregated drugs pr. target
#' @param indication_label column name for drug indications
#'
#'
#' @export
targeted_drugs_summarise <- function(
  candidate_drugs = NULL,
  link_label = "DRUGS_ON_LABEL",
  indication_label = "DRUGS_ON_LABEL_INDICATIONS"){

  invisible(assertthat::assert_that(!is.null(candidate_drugs)))
  invisible(assertthat::assert_that(is.data.frame(candidate_drugs)))

  assertable::assert_colnames(candidate_drugs,
                              c('symbol', 'nci_concept_display_name',
                              'drug_link', 'max_phase',
                              'drug_indication_label'),
                              only_colnames = F, quiet = T)

  candidate_drugs <- as.data.frame(
    candidate_drugs %>%
      dplyr::select(.data$symbol, .data$nci_concept_display_name,
                    .data$drug_link, .data$max_phase, .data$drug_indication_label) %>%
      tidyr::separate_rows(.data$drug_indication_label, sep = "\\|") %>%
      dplyr::arrange(.data$symbol, dplyr::desc(.data$max_phase)) %>%
      dplyr::distinct() %>%
      dplyr::group_by(.data$symbol, .data$nci_concept_display_name, .data$drug_link) %>%
      dplyr::summarise(drug_indication_label =
                         paste(unique(.data$drug_indication_label),
                               collapse = ", "),
                       .groups = "drop") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(drug_indication_label =
                      paste0(.data$drug_indication_label,
                             " (<b>", .data$nci_concept_display_name, "</b>)")) %>%
      dplyr::group_by(.data$symbol) %>%
      dplyr::summarise(!!rlang::sym(link_label) := paste0(.data$drug_link, collapse = ", "),
                       !!rlang::sym(indication_label) := paste0(.data$drug_indication_label,
                                                         collapse = "; "),
                       .groups = "drop")
  )

  return(candidate_drugs)


}

log4r_info <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  log4r::info(log4r_logger, msg)
}

log4r_debug <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  log4r::debug(log4r_logger, msg)
}

log4r_warn <- function(msg) {
  log4r_logger <- getOption("PCGRR_LOG4R_LOGGER")
  log4r::warn(log4r_logger, msg)
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
      "Please install with:\n'BiocManager::install(\"{pkg}\")'"
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
