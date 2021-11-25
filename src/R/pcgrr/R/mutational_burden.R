#' Function that tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#'
#' @return pcg_report_data data frame with all report elements
#' @export
generate_report_data_tmb <- function(sample_calls,
                                     pcgr_data,
                                     sample_name, pcgr_config) {

  tmb_consequence_pattern <-
    "^(stop_|start_lost|frameshift_|missense_|synonymous_|inframe_)"

  log4r_info("------")
  log4r_info(paste0("Calculating tumor mutational burden"))

  pcg_report_tmb <-
    pcgrr::init_report(config = pcgr_config,
                       class = "tmb",
                       pcgr_data = pcgr_data)

  if (pcgr_config[["tmb"]][["algorithm"]] == "nonsyn") {
    tmb_consequence_pattern <- "^(missense_)"
    pcg_report_tmb[["tmb"]][["algorithm"]] <-
      pcgr_config[["tmb"]][["algorithm"]]
  }

  pcg_report_tmb[["eval"]] <- TRUE

  if(NROW(sample_calls) > 0){
    pcg_report_tmb[["v_stat"]][["n_tmb"]] <-
      sample_calls %>%
      dplyr::filter(
        stringr::str_detect(.data$CONSEQUENCE, tmb_consequence_pattern)) %>%
      nrow()
  }

  if (pcg_report_tmb[["v_stat"]][["n_tmb"]] > 0 &
      pcgr_config$assay_props$target_size_mb > 0) {
    pcg_report_tmb[["v_stat"]][["tmb_estimate"]] <-
      round(
        as.numeric(pcg_report_tmb[["v_stat"]][["n_tmb"]] /
                     pcg_report_tmb[["v_stat"]][["target_size_mb"]]),
        digits = 2)
  }
  log4r_info(
    paste0("Number of variants for mutational burden analysis: ",
           pcg_report_tmb[["v_stat"]][["n_tmb"]]))
  log4r_info(
    paste0("Coding target size sequencing assay (Mb): ",
           pcg_report_tmb[["v_stat"]][["target_size_mb"]]))
  log4r_info(
    paste0("Estimated mutational burden: ",
           pcg_report_tmb[["v_stat"]][["tmb_estimate"]],
           " mutations/Mb"))

  return(pcg_report_tmb)
}

#' Function that makes a plot with TMB boxplots for TCGA cohorts, highlighting
#' the TMB estimate for a given sample and the cohort/primary site of interest
#'
#' @param tcga_tmb data frame with TMB estimates for TCGA samples
#' @param p_site primary tumor_site (sample)
#' @param tmb_estimate estimate of mutational burden (sample)
#' @param algorithm type of TMB algorithm
#'
#'
#' @export
plot_tmb_primary_site_tcga <- function(tcga_tmb, p_site = "Liver",
                                       tmb_estimate = 5,
                                       algorithm = "all_coding") {


  tmb_site_colors <- data.frame(primary_site =
                                  unique(tcga_tmb$primary_site),
                                stringsAsFactors = F) %>%
    dplyr::filter(!is.na(.data$primary_site))
  tmb_site_colors$color <- "#f0f0f0"
  tmb_site_colors <-
    dplyr::mutate(tmb_site_colors,
                  color = dplyr::if_else(
                    .data$primary_site == p_site,
                    pcgrr::color_palette[["tier"]][["values"]][1], .data$color))

  tmb_site_color_vec <- tmb_site_colors$color
  names(tmb_site_color_vec) <- tmb_site_colors$primary_site

  tcga_tmb <- tcga_tmb %>%
    dplyr::filter(!is.na(.data$primary_site)) %>%
    dplyr::select(.data$primary_site, .data$tmb_log10, .data$tmb)

  if (algorithm == "nonsyn") {
    tcga_tmb <- tcga_tmb %>%
      dplyr::filter(!is.na(.data$primary_site)) %>%
      dplyr::select(.data$primary_site, .data$tmb_ns_log10, .data$tmb_ns) %>%
      dplyr::rename(tmb = .data$tmb_ns, tmb_log10 = .data$tmb_ns_log10)
  }

  tmb_plot_site <-
    ggplot2::ggplot(data = tcga_tmb) +
    ggplot2::geom_boxplot(mapping =
                            ggplot2::aes(x = stats::reorder(.data$primary_site, .data$tmb_log10,
                                                     FUN = stats::median),
                                         y = .data$tmb, fill = .data$primary_site)) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(trans = scales::log_trans(base = 10),
                                breaks = c(0.01, 1, 10, 100, 1000),
                                labels = c("0.01", "1", "10", "100", "1000")) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = tmb_site_color_vec) +
    ggplot2::xlab("Primary Site") +
    ggplot2::ylab("Tumor mutational burden (mutations/mb)") +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "none",
                   axis.text.x =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 14),
                   axis.title.x =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 16, vjust = -1.5),
                   axis.text.y =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 14),
                   plot.margin =
                     (grid::unit(c(0.5, 2, 2, 0.5), "cm")),
                   legend.text = ggplot2::element_text(family = "Helvetica",
                                                       size = 14)) +
    ggplot2::geom_hline(yintercept = as.numeric(tmb_estimate), size = 0.9,
                        linetype = 4,
                        colour = pcgrr::color_palette[["tier"]][["values"]][1])


  return(tmb_plot_site)

}
