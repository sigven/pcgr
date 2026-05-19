#' Function that assigns one of six mutation types to a list of mutations
#'
#' @param var_df data frame with variants
#'
#' @return var_df
#'
#' @export
assign_mutation_type <- function(var_df) {

  invisible(
    assertthat::assert_that(
      is.data.frame(var_df),
      msg = "Argument 'var_df' must be a valid data.frame"))
  assertable::assert_colnames(var_df, c("VARIANT_CLASS", "REF", "ALT"),
                              only_colnames = F, quiet = T)
  var_df <- var_df |>
    dplyr::mutate(
      MUTATION_TYPE =
        dplyr::case_when(VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "A" ~ "C>T:G>A",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "T" ~ "C>T:G>A",
                         VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "C" ~ "C>G:G>C",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "G" ~ "C>G:G>C",
                         VARIANT_CLASS == "SNV" & REF == "C" &
                           ALT == "A" ~ "C>A:G>T",
                         VARIANT_CLASS == "SNV" & REF == "G" &
                           ALT == "T" ~ "C>A:G>T",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "G" ~ "A>G:T>C",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "C" ~ "A>G:T>C",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "T" ~ "A>T:T>A",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "A" ~ "A>T:T>A",
                         VARIANT_CLASS == "SNV" & REF == "A" &
                           ALT == "C" ~ "A>C:T>G",
                         VARIANT_CLASS == "SNV" & REF == "T" &
                           ALT == "G" ~ "A>C:T>G",
                         FALSE ~ as.character(NA)))
  return(var_df)
}

#' Function that generates a VAF distribution plot for a given PCGR report object
#'
#' @param report PCGR report object
#' @param font_size font size for plot text elements
#' @param font_family font family for plot text elements
#' @return vaf_plot_plotly plotly object with VAF distribution plot
#'
#' @export
#'
#'
vaf_plot <- function(
    report = NULL,
    font_size = 12,
    font_family = "Helvetica") {

  invisible(assertthat::assert_that(
    !is.null(report),
    msg = "Argument 'report' must be provided"))
  invisible(assertthat::assert_that(
    is.list(report),
    msg = "Argument 'report' must be a valid PCGR report object"))
  invisible(assertthat::assert_that(
    !is.null(report$content$snv_indel$callset$variant),
    msg = "Argument 'report' must contain a valid callset"))

  vaf_dist_tumor <- af_distribution(
    var_df = report$content$snv_indel$callset$variant)

  assertable::assert_colnames(
    vaf_dist_tumor,
    c("bin_start","bin_end","Count","VARIANT_CLASS"),
    only_colnames = F,
    quiet = T
  )

  vaf_dist_tumor$VARIANT_CLASS <- factor(
    vaf_dist_tumor$VARIANT_CLASS,
    levels = c("SNV","deletion","insertion","indel","substitution"))

  vaf_plot <- ggplot2::ggplot(data = vaf_dist_tumor) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(
        x = .data$bin_start,
        y = .data$Count,
        fill = .data$VARIANT_CLASS),
      stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::ylab("Number of variants") +
    ggplot2::xlab("Variant allelic fraction - tumor") +
    ggplot2::scale_fill_manual(
      values = color_palette$multi$values) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      #legend.position = "",
      #axis.text.x = element_blank(),
      axis.text.x = ggplot2::element_text(
        family = font_family,
        size = font_size,
        vjust = -0.1),
      axis.title.x = ggplot2::element_text(
        family = font_family,
        size = font_size,
        vjust = -0.5),
      axis.text.y = ggplot2::element_text(
        family = "Helvetica",
        size = font_size),
      axis.title.y = ggplot2::element_text(
        family = "Helvetica",
        size = font_size,
        vjust = 1.5),
      plot.margin = (grid::unit(c(0.5, 0.5, 1.5, 0.5), "cm")),
      legend.text = ggplot2::element_text(
        family =
          font_family,
        size = font_size))

  vaf_plot_plotly <- plotly::ggplotly(vaf_plot)
  vaf_plot_plotly$x$layout$legend$title$text <- ""
  vaf_plot_plotly <- vaf_plot_plotly |>
    plotly::layout(
      legend = list(orientation = "h",
                    x = 0.5, xanchor = "center",
                    y = -0.25, yanchor = "top"),
      margin = list(b = 110, t = 30))

  return(vaf_plot_plotly)

}

#' Clean gnomAD VCF annotations
#'
#' Function that cleans and processes gnomAD annotations in the
#' raw VCF-annotated callset for SNVs/InDels
#'
#' @param var_df data frame with variants
#' @return var_df with cleaned gnomAD annotations
#'
#' @export
#'
clean_gnomad_annotations <- function(
    var_df = NULL) {

  if (is.null(var_df) | !is.data.frame(var_df)) {
    return(var_df)
  }

  if (NROW(var_df) == 0) {
    return(var_df)
  }

  ## gnomAD full allele frequency (v4.1) - predisposition targets
  if ("gALL_GRPMAX" %in% colnames(var_df)) {

    new_cols <-
      c("gnomADj_AC_GRPMAX",
        "gnomADj_AN_GRPMAX",
        "gnomADj_NHOMALT_GRPMAX",
        "gnomADj_POP_GRPMAX",
        "gnomADj_FAF95_GRPMAX",
        "gnomADj_FAF95_POP_GRPMAX")

    gnomad_full_cols <-
      c("tmp_gALL_AC_joint",
        "tmp_gALL_AN_joint",
        "tmp_gALL_AC_GRPMAX_joint",
        "tmp_gALL_AN_GRPMAX_joint",
        "tmp_gALL_NHOMALT_GRPMAX_joint",
        "tmp_gALL_POP_GRPMAX_joint",
        "tmp_gALL_FAF95_GRPMAX_joint",
        "tmp_gALL_FAF95_POP_GRPMAX_joint",
        "tmp_gALL_AC_exomes",
        "tmp_gALL_AN_exomes",
        "tmp_gALL_AC_GRPMAX_exomes",
        "tmp_gALL_AN_GRPMAX_exomes",
        "tmp_gALL_NHOMALT_GRPMAX_exomes",
        "tmp_gALL_POP_GRPMAX_exomes",
        "tmp_gALL_FAF95_GRPMAX_exomes",
        "tmp_gALL_FAF95_POP_GRPMAX_exomes",
        "tmp_gALL_AC_genomes",
        "tmp_gALL_AN_genomes",
        "tmp_gALL_AC_GRPMAX_genomes",
        "tmp_gALL_AN_GRPMAX_genomes",
        "tmp_gALL_NHOMALT_GRPMAX_genomes",
        "tmp_gALL_POP_GRPMAX_genomes",
        "tmp_gALL_FAF95_GRPMAX_genomes",
        "tmp_gALL_FAF95_POP_GRPMAX_genomes",
        "tmp_gALL_not_called_in_genomes",
        "tmp_gALL_not_called_in_exomes",
        "tmp_gALL_exomes_filters",
        "tmp_gALL_genomes_filters")

    for (col in gnomad_full_cols) {
      if (col %in% colnames(var_df)) {
        var_df[[col]] <- NULL
      }
    }

    for (col in new_cols) {
      if (col %in% colnames(var_df)) {
        var_df[[col]] <- NULL
      }
    }

    var_df <-
      var_df |>
      tidyr::separate(
        col = "gALL_GRPMAX",
        into = gnomad_full_cols,
        sep = "\\|",
        remove = TRUE,
        fill = "right",
        extra = "drop"
      ) |>
      dplyr::mutate(
        gnomADj_AC_GRPMAX =
          as.integer(.data$tmp_gALL_AC_GRPMAX_joint),
        gnomADj_AN_GRPMAX =
          as.integer(.data$tmp_gALL_AN_GRPMAX_joint),
        gnomADj_NHOMALT_GRPMAX =
          as.integer(.data$tmp_gALL_NHOMALT_GRPMAX_joint),
        gnomADj_POP_GRPMAX = as.character(
          .data$tmp_gALL_POP_GRPMAX_joint),
        gnomADj_FAF_GRPMAX = as.numeric(
          .data$tmp_gALL_FAF95_GRPMAX_joint),
        gnomADj_FAF_POP_GRPMAX = as.character(
          .data$tmp_gALL_FAF95_POP_GRPMAX_joint)
      )

    for (c in gnomad_full_cols) {
      if (c %in% colnames(var_df)) {
        var_df[[c]] <- NULL
      }
    }
  }

  # check if gnomAD non-cancer filter allele frequency is available
  if ("gNC_FAF" %in% colnames(var_df)) {
    for (pop in c("GRPMAX","AFR","AMR","NFE","EAS","SAS","GLOBAL")) {
      tag <- paste0("tmp_gNC_FAF_", pop)
      if (tag %in% colnames(var_df)) {
        var_df[[tag]] <- NULL
      }
    }
    var_df <-
      var_df |>
      tidyr::separate(
        col = "gNC_FAF",
        into = c("tmp_gNC_FAF_GRPMAX",
                 "tmp_gNC_FAF_GLOBAL",
                 "tmp_gNC_FAF_AFR",
                 "tmp_gNC_FAF_AMR",
                 "tmp_gNC_FAF_EAS",
                 "tmp_gNC_FAF_FIN",
                 "tmp_gNC_FAF_NFE",
                 "tmp_gNC_FAF_SAS",
                 "tmp_gNC_VAR_FILTER"),
        sep = "\\|",
        remove = T,
        fill = "right",
        extra = "drop"
      )

    ## convert to numeric; non-numeric tokens (e.g. ".", "") become NA intentionally
    for (pop in c("GRPMAX","AFR","AMR","EAS","FIN","NFE","SAS","GLOBAL")) {
      tag <- paste0("tmp_gNC_FAF_", pop)
      if (tag %in% colnames(var_df)) {
        var_df <-
          var_df |>
          dplyr::mutate(
            !!rlang::sym(tag) := suppressWarnings(as.numeric(!!rlang::sym(tag)))
          )
      }
    }

    var_df <-
      var_df |>
      grpmax_faf_nc_gnomad()

    var_df[['gnomAD_NC_VAR_FILTER']] <-
      var_df[['tmp_gNC_VAR_FILTER']]
    var_df[['tmp_gNC_VAR_FILTER']] <- NULL
    for (pop in c("AFR","AMR","EAS","FIN","NFE","SAS","GLOBAL","GRPMAX")) {
      tag <- paste0("tmp_gNC_FAF_", pop)
      if (tag %in% colnames(var_df)) {
        var_df[[tag]] <- NULL
      }
    }

  }


  ## check if gnomAD non-cancer allele count/allele number is available
  ## - make popmax and global allele frequency variables
  ## - only consider populations with AN >= 4000

  var_df[["gnomAD_NC_AF_POPMAX"]] <- -1000
  var_df[["gnomAD_NC_AC_POPMAX"]] <- -1000
  var_df[["gnomAD_NC_POP_POPMAX"]] <- NA
  var_df[["gnomAD_NC_AF_GLOBAL"]] <- -1000
  var_df[["gnomAD_NC_AN_GLOBAL"]] <- -1000

  for (pop in c('AMR','AFR','EAS','FIN','NFE','SAS','')) {
    gnc_tag <- paste0("gNC_", pop)
    if (pop == '') {
      gnc_tag <- "gNC"
    }

    if (gnc_tag %in% colnames(var_df)) {

      tags <- list()
      for (t in c('AC','AN','NHOMALT','AF')) {
        tags[[t]] <- gsub('__','_',paste0("tmp_gnomAD_non_cancer_", pop, "_", t))
      }

      ## if the gnomAD tag is present, split it into
      ## the expanded tag set
      var_df <- var_df |>
        tidyr::separate(
          !!rlang::sym(gnc_tag),
          into = c(tags[['AC']],
                   tags[['AN']],
                   tags[['NHOMALT']]),
          sep = "\\|",
          remove = TRUE,
          fill = "right",
          extra = "drop"
        ) |>
        dplyr::mutate(
          !!rlang::sym(tags[['AN']]) := suppressWarnings(
            as.integer(!!rlang::sym(tags[['AN']]))
          ),
          !!rlang::sym(tags[['AC']]) := suppressWarnings(
            as.integer(!!rlang::sym(tags[['AC']]))
          ),
          !!rlang::sym(tags[['NHOMALT']]) := suppressWarnings(
            as.integer(!!rlang::sym(tags[['NHOMALT']]))
          ),
          !!rlang::sym(tags[['AF']]) := as.numeric(
            as.numeric(!!rlang::sym(tags[['AC']])) /
              as.integer(!!rlang::sym(tags[['AN']])))

        )

      pop2 <- pop
      if (pop == "") {
        pop2 <- "GLOBAL"
        var_df[["gnomAD_NC_AF_GLOBAL"]] <-
          var_df[[tags[['AF']]]]
        var_df[["gnomAD_NC_AN_GLOBAL"]] <-
          var_df[[tags[['AN']]]]
      }
      var_df[!is.na(var_df[[tags[['AF']]]]) &
               !is.na(var_df[[tags[['AN']]]]) &
               var_df[[tags[['AF']]]] >=
               var_df$gnomAD_NC_AF_POPMAX &
               var_df[[tags[['AN']]]] >= 4000, "gnomAD_NC_POP_POPMAX"] <-
        tolower(as.character(pop2))

      var_df[!is.na(var_df[[tags[['AF']]]]) &
               !is.na(var_df[[tags[['AN']]]]) &
               var_df[[tags[['AF']]]] >=
               var_df$gnomAD_NC_AF_POPMAX &
               var_df[[tags[['AN']]]] >= 4000, "gnomAD_NC_AF_POPMAX"] <- as.numeric(
                 var_df[!is.na(var_df[[tags[['AF']]]]) &
                          !is.na(var_df[[tags[['AN']]]]) &
                          var_df[[tags[['AF']]]] >=
                          var_df$gnomAD_NC_AF_POPMAX &
                          var_df[[tags[['AN']]]] >= 4000, tags[['AF']]])

      var_df[!is.na(var_df[[tags[['AF']]]]) &
               !is.na(var_df[[tags[['AN']]]]) &
               var_df[[tags[['AF']]]] >=
               var_df$gnomAD_NC_AF_POPMAX &
               var_df[[tags[['AN']]]] >= 4000, "gnomAD_NC_AC_POPMAX"] <- as.integer(
                 var_df[!is.na(var_df[[tags[['AF']]]]) &
                          !is.na(var_df[[tags[['AN']]]]) &
                          var_df[[tags[['AF']]]] >=
                          var_df$gnomAD_NC_AF_POPMAX &
                          var_df[[tags[['AN']]]] >= 4000, tags[['AC']]])

      var_df[[tags[['AN']]]] <- NULL
      var_df[[tags[['AC']]]] <- NULL
      var_df[[tags[['NHOMALT']]]] <- NULL
      var_df[[tags[['AF']]]] <- NULL

    }
  }
  ## convert -1000 to NA for gnomAD NC AF/AC/AN variables
  if ("gnomAD_NC_AC_POPMAX" %in% colnames(var_df) &
     "gnomAD_NC_AF_POPMAX" %in% colnames(var_df) &
     "gnomAD_NC_AN_GLOBAL" %in% colnames(var_df) &
     "gnomAD_NC_AF_GLOBAL" %in% colnames(var_df)) {
    var_df <- var_df |>
      dplyr::mutate(
        gnomAD_NC_AF_POPMAX = dplyr::case_when(
          .data$gnomAD_NC_AF_POPMAX == -1000 ~ as.numeric(NA),
          TRUE ~ .data$gnomAD_NC_AF_POPMAX
        )) |>
      dplyr::mutate(
        gnomAD_NC_AC_POPMAX = dplyr::case_when(
          .data$gnomAD_NC_AC_POPMAX == -1000 ~ as.integer(NA),
          TRUE ~ .data$gnomAD_NC_AC_POPMAX
        )) |>
      dplyr::mutate(
        gnomAD_NC_AN_GLOBAL = dplyr::case_when(
          .data$gnomAD_NC_AN_GLOBAL == -1000 ~ as.integer(NA),
          TRUE ~ .data$gnomAD_NC_AN_GLOBAL
        )
      ) |>
      dplyr::mutate(
        gnomAD_NC_AF_GLOBAL = dplyr::case_when(
          .data$gnomAD_NC_AF_GLOBAL == -1000 ~ as.numeric(NA),
          TRUE ~ .data$gnomAD_NC_AF_GLOBAL
        )
      )
  }

  return(var_df)

}
