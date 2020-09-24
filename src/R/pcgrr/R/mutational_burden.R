#' Function that tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#'
#' @return pcg_report_data data frame with all report elements
#'
generate_report_data_tmb <- function(sample_calls, pcgr_data, sample_name, pcgr_config) {

  tmb_consequence_pattern <- "^(stop_|start_lost|frameshift_|missense_|synonymous_|inframe_)"
  rlogging::message("------")
  rlogging::message(paste0("Calculating tumor mutational burden"))

  pcg_report_tmb <- pcgrr::init_report(pcgr_config, sample_name = sample_name,
                                           class = "tmb", pcgr_data = pcgr_data)

  pcg_report_tmb[["eval"]] <- TRUE
  pcg_report_tmb[["variant_statistic"]][["n_tmb"]] <- sample_calls %>%
    dplyr::filter(stringr::str_detect(CONSEQUENCE, tmb_consequence_pattern)) %>% nrow()
  if (pcg_report_tmb[["variant_statistic"]][["n_tmb"]] > 0 & pcgr_config$assay_properties$target_size_mb > 0) {
    pcg_report_tmb[["variant_statistic"]][["tmb_estimate"]] <-
      round(as.numeric(pcg_report_tmb[["variant_statistic"]][["n_tmb"]] /
                         pcg_report_tmb[["variant_statistic"]][["target_size_mb"]]), digits = 2)
    # if (pcg_report_tmb[["variant_statistic"]][["tmb_estimate"]] <= pcgr_config$tmb$tmb_low_limit) {
    #   pcg_report_tmb[["variant_statistic"]][["tmb_tertile"]] <-
    #     paste0("TMB - Low\n(0 - ", pcgr_config$tmb$tmb_low_limit, " mutations/Mb)")
    # }
    # if (pcg_report_tmb[["variant_statistic"]][["tmb_estimate"]] > pcgr_config$tmb$tmb_low_limit &
    #     pcg_report_tmb[["variant_statistic"]][["tmb_estimate"]] <= pcgr_config$tmb$tmb_intermediate_limit) {
    #   pcg_report_tmb[["variant_statistic"]][["tmb_tertile"]] <-
    #     paste0("TMB - Intermediate\n(", pcgr_config$tmb$tmb_low_limit, " - ",
    #            pcgr_config$tmb$tmb_intermediate_limit, " mutations/Mb)")
    # }
    # if (pcg_report_tmb[["variant_statistic"]][["tmb_estimate"]] > pcgr_config$tmb$tmb_intermediate_limit) {
    #   pcg_report_tmb[["variant_statistic"]][["tmb_tertile"]] <-
    #     paste0("TMB - High\n( > ", pcgr_config$tmb$tmb_intermediate_limit, " mutations/Mb)")
    # }
  }
  rlogging::message(paste0("Number of variants for mutational burden analysis: ", pcg_report_tmb[["variant_statistic"]][["n_tmb"]]))
  rlogging::message(paste0("Coding target size sequencing assay (Mb): ", pcg_report_tmb[["variant_statistic"]][["target_size_mb"]]))
  rlogging::message(paste0("Estimated mutational burden: ", pcg_report_tmb[["variant_statistic"]][["tmb_estimate"]], " mutations/Mb"))
  #rlogging::message(paste0("Mutational burden tertile: ", stringr::str_replace(pcg_report_tmb[["variant_statistic"]][["tmb_tertile"]], "\n", " ")))

  return(pcg_report_tmb)
}

plot_tmb_primary_site_tcga <- function(tcga_tmb, p_site = "Liver", tmb_estimate = 5, tmb_high = 20) {


  tmb_site_colors <- data.frame(primary_site = unique(tcga_tmb$primary_site), stringsAsFactors = F) %>%
    dplyr::filter(!is.na(primary_site))
  tmb_site_colors$color <- "#f0f0f0"
  tmb_site_colors <- dplyr::mutate(tmb_site_colors,
                                   color = dplyr::if_else(primary_site == p_site,
                                                          pcgrr::color_palette[['tier']][['values']][1], color))

  tmb_site_color_vec <- tmb_site_colors$color
  names(tmb_site_color_vec) <- tmb_site_colors$primary_site

  tcga_tmb <- tcga_tmb %>% dplyr::filter(!is.na(primary_site))

  tmb_plot_site <-
    ggplot2::ggplot(data = tcga_tmb) +
    ggplot2::geom_boxplot(mapping = ggplot2::aes(x = reorder(primary_site, tmb_log10, FUN = median),
                                                 y = tmb, fill = primary_site)) +
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
                   axis.text.x = ggplot2::element_text(family = "Helvetica", size = 14),
                   axis.title.x = ggplot2::element_text(family = "Helvetica", size = 16, vjust = -1.5),
                   axis.text.y = ggplot2::element_text(family = "Helvetica", size = 14),
                   plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm")),
                   legend.text = ggplot2::element_text(family = "Helvetica", size = 14)) +
    ggplot2::geom_hline(yintercept = as.numeric(tmb_estimate), size = 0.9,
                        linetype = 4, colour = pcgrr::color_palette[['tier']][['values']][1])
    #ggplot2::geom_rect(ggplot2::aes(ymin = tmb_high, ymax = 1020, xmin = -Inf, xmax = Inf),
                       #fill = "gray", alpha = 0.01)

  return(tmb_plot_site)

}
