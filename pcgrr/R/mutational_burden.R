#' Function that reads TSV file with TMB estimates from sample
#'
#' @param settings PCGR run/configuration settings
#'
#' @return tmb_estimate
#' @export
generate_report_data_tmb <- function(settings = NULL) {

  tmb_data <- data.frame()
  if(settings[['molecular_data']][['fname_tmb_tsv']] != "None" &
     file.exists(
       settings[['molecular_data']][['fname_tmb_tsv']]
     )){

    tmb_data <- readr::read_tsv(
      settings[['molecular_data']][['fname_tmb_tsv']],
      show_col_types = F
    )

  }

  tmb_rep <- list()
  tmb_rep[["eval"]] <- TRUE
  tmb_rep[["sample_estimate"]] <- tmb_data

  return(tmb_rep)
}

#' Function that makes a plot with TMB boxplots for reference cohorts, highlighting
#' the TMB estimate for a given sample and the cohort/primary site of interest
#'
#' @param tmb_reference data frame with TMB estimates for reference samples (e.g. TCGA)
#' @param p_site primary tumor_site (sample)
#' @param tmb_estimates data frame with estimates of mutational burden (sample)
#' @param tmb_display_type TMB estimation used for display (e.g. missense_only)
#' @param tumor_only logical, if TRUE color using tumor-only color
#'
#'
#' @export
plot_tmb_primary_site_tcga <- function(
    tmb_reference,
    p_site = "Liver",
    tmb_estimates = NULL,
    tmb_display_type = "missense_only",
    tumor_only = FALSE) {


  report_color <- pcgrr::color_palette$report_color$values[1]
  if(tumor_only == TRUE){
    report_color <- pcgrr::color_palette$report_color$values[2]
  }

  assertable::assert_colnames(
    tmb_reference,
    c("PRIMARY_SITE",
      "TMB_CODING_AND_SILENT",
      "TMB_CODING_NON_SILENT",
      "TMB_MISSENSE_ONLY",
      "PRIMARY_SITE",
      "PRIMARY_DIAGNOSIS_VERY_SIMPLIFIED"),
    only_colnames = F,
    quiet = T
  )

  tmb_estimate <- tmb_estimates |>
    dplyr::filter(
      .data$tmb_measure == "TMB_missense_only")

  if(tmb_display_type == "coding_and_silent"){
    tmb_estimate <- tmb_estimates |>
      dplyr::filter(
        .data$tmb_measure == "TMB_coding_and_silent")
  }
  if(tmb_display_type == "coding_non_silent"){
    tmb_estimate <- tmb_estimates |>
      dplyr::filter(
        .data$tmb_measure == "TMB_coding_non_silent")
  }


  tmb_site_colors <- data.frame(
    PRIMARY_SITE =
      unique(tmb_reference$PRIMARY_SITE),
    stringsAsFactors = F) |>
    dplyr::filter(!is.na(.data$PRIMARY_SITE))
  tmb_site_colors$color <- "#BABABA"
  tmb_site_colors <-
    dplyr::mutate(
      tmb_site_colors,
      color = dplyr::if_else(
        .data$PRIMARY_SITE == p_site,
        report_color,
        .data$color))

  tmb_site_color_vec <- tmb_site_colors$color
  names(tmb_site_color_vec) <- tmb_site_colors$PRIMARY_SITE

  tmb_reference <- tmb_reference |>
    dplyr::filter(!is.na(.data$PRIMARY_SITE)) |>
    dplyr::rename(
      SUBTYPE = "PRIMARY_DIAGNOSIS_VERY_SIMPLIFIED") |>
    dplyr::mutate(
      TMB = dplyr::case_when(
        tmb_display_type == "coding_and_silent" ~
          as.numeric(.data$TMB_CODING_AND_SILENT) + 0.001,
        tmb_display_type == "coding_non_ssilent" ~
          as.numeric(.data$TMB_CODING_NON_SILENT) + 0.001,
        tmb_display_type == "missense_only" ~
          as.numeric(.data$TMB_MISSENSE_ONLY) + 0.001,
        TRUE ~ as.numeric(NA)
      )) |>
    dplyr::select(
      c("PRIMARY_SITE", "TMB",
        "SUBTYPE"))

  median_tmb <- as.data.frame(tmb_reference |>
    dplyr::group_by(.data$PRIMARY_SITE) |>
    dplyr::summarize(
      TMB_MEDIAN = stats::median(.data$TMB, na.rm = T),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$TMB_MEDIAN)))

  tmb_reference$PRIMARY_SITE <- factor(
    tmb_reference$PRIMARY_SITE,
    levels = median_tmb$PRIMARY_SITE)

  tmb_plot_site <-
    ggplot2::ggplot(
      data = tmb_reference, mapping =
        ggplot2::aes(x = .data$PRIMARY_SITE,
                     y = .data$TMB,
                     fill = .data$PRIMARY_SITE)) +
    ggplot2::geom_boxplot(width = 0.4) +
    ggplot2::geom_violin() +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(
      trans = scales::log_trans(base = 10),
      breaks = c(0.01, 0.1, 1, 10, 100, 1000),
      labels = c("-2", "-1", "0", "1", "2", "3")) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = tmb_site_color_vec) +
    #ggplot2::xlab("Primary Site") +
    ggplot2::ylab("log10 tumor mutational burden (mutations/mb)") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x =
        ggplot2::element_text(family = "Helvetica",
                              size = 12),
      axis.title.x =
        ggplot2::element_text(family = "Helvetica",
                              size = 12, vjust = -1.5),
      axis.text.y =
        ggplot2::element_text(family = "Helvetica",
                              size = 12),
      plot.margin =
        (grid::unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
      legend.text = ggplot2::element_blank())

  if (NROW(tmb_estimate) == 1) {
    tmb_plot_site <- tmb_plot_site +
      ggplot2::geom_hline(
        yintercept = as.numeric(tmb_estimate$tmb_estimate),
        linewidth = 0.6,
        linetype = "dashed",
        colour = report_color)
  }


  return(tmb_plot_site)

}
