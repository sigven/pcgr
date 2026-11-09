## PDF report (Quarto/Typst) - static (ggplot2) figures and plot settings

#' Font settings for (static) ggplot2 figures of the PDF report - generic
#' sans-serif font family and smaller text than the HTML report plots
#' (which set 'Helvetica', not available on all systems)
#'
#' @param p ggplot2 object
#' @param base_size base text size (pt)
#'
#' @return ggplot2 object
#' @keywords internal
pdf_plot_theme <- function(p, base_size = 9) {
  p + ggplot2::theme(
    text = ggplot2::element_text(family = "sans", size = base_size),
    axis.text.x = ggplot2::element_text(family = "sans", size = base_size,
                                        face = "plain"),
    axis.text.y = ggplot2::element_text(family = "sans", size = base_size),
    axis.title.x = ggplot2::element_text(family = "sans", size = base_size),
    legend.text = ggplot2::element_text(family = "sans", size = base_size),
    plot.margin = ggplot2::margin(4, 6, 4, 4))
}

#' Donut chart of variant statistics (static, ggplot2) for the PDF report -
#' counterpart of plotly_pie_chart() in the HTML report
#'
#' @param df_variant_stats data frame with 'N', 'Pct' and the category
#' column, from stats_type_snv_indel()
#' @param category category column (e.g. 'CODING_STATUS')
#' @param title chart title
#' @param center_label text in the center of the donut (default: total N)
#' @param show_counts show counts (N) with the legend labels
#' @param color_palette PCGR color palette
#'
#' @return ggplot2 object
#' @keywords internal
pdf_plot_donut <- function(df_variant_stats,
                          category = "CODING_STATUS",
                          title = NULL,
                          center_label = NULL,
                          show_counts = TRUE,
                          color_palette = pcgrr::color_palette) {
  d <- df_variant_stats
  d$LABEL <- stringr::str_replace_all(as.character(d[[category]]), "_", " ")
  d <- d[order(-d$N), , drop = FALSE]
  ## legend labels truncated, so that they fit the fixed legend margin
  d$LABEL <- stringr::str_trunc(d$LABEL, 30)
  if (show_counts) d$LABEL <- paste0(d$LABEL, " (", d$N, ")")
  d$LABEL <- factor(d$LABEL, levels = d$LABEL)
  if (is.null(center_label)) center_label <- paste0("N = ", sum(d$N))
  d$ymax <- cumsum(d$N)
  d$ymin <- d$ymax - d$N
  d$ymid <- (d$ymax + d$ymin) / 2
  d$PCT_LABEL <- ifelse(d$Pct >= 6, paste0(round(d$Pct), "%"), "")
  cols <- utils::head(color_palette$multi$values, NROW(d))

  ggplot2::ggplot(d) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = 2, xmax = 3.4, ymin = .data$ymin,
                   ymax = .data$ymax, fill = .data$LABEL),
      colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(x = 2.7, y = .data$ymid, label = .data$PCT_LABEL),
      colour = "white", fontface = "bold", size = 2.7) +
    ggplot2::annotate(
      "text", x = 0.6, y = 0, label = center_label,
      fontface = "bold", size = 3.2, colour = color_palette$bg_dark) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(0.6, 3.4)) +
    ggplot2::scale_fill_manual(values = cols) +
    ## a title line is always reserved, so that all donuts (with the same
    ## figure size) are drawn at the same size
    ggplot2::labs(title = if (is.null(title)) " " else title) +
    ggplot2::theme_void(base_size = 7) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 8, colour = color_palette$bg_dark,
        margin = ggplot2::margin(b = 2)),
      ## legend placed in a fixed right margin (not next to the panel), so
      ## that the donut size does not depend on the legend width
      legend.position = "inside",
      legend.position.inside = c(1.04, 0.5),
      legend.justification.inside = c(0, 0.5),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 6.5),
      legend.key.size = grid::unit(0.28, "cm"),
      plot.margin = ggplot2::margin(2, 125, 2, 2))
}

#' Donut charts of SNV/InDel variant statistics (coding status, variant
#' type, consequence types) for the PDF report
#'
#' @param report PCGR report object
#'
#' @return list of ggplot2 objects (coding consequence chart only when
#' coding variants are present)
#' @keywords internal
pdf_plot_varstats <- function(report) {
  variant_stats <- pcgrr::stats_type_snv_indel(
    var_df = report$content$snv_indel$callset$variant_display,
    pct_other_limit = 3.5)
  plots <- list(
    coding = pdf_plot_donut(
      variant_stats[["coding"]], "CODING_STATUS", "Coding effect"),
    variant_class = pdf_plot_donut(
      variant_stats[["variant_class"]], "VARIANT_CLASS", "Variant type"),
    consequence = pdf_plot_donut(
      variant_stats[["consequence"]], "CONSEQUENCE",
      "Consequence type - all variants"))
  if (NROW(variant_stats[["consequence_coding"]]) > 0) {
    plots[["consequence_coding"]] <- pdf_plot_donut(
      variant_stats[["consequence_coding"]], "CONSEQUENCE",
      "Consequence type - coding variants")
  }
  plots
}

#' Copy number plot settings (tumor ploidy and effective amplification /
#' gain / deletion thresholds) for the PDF report - as used for the copy
#' number profile plots of the HTML report
#'
#' @param report PCGR report object
#'
#' @return list with tumor_ploidy, amp_threshold_effective, threshold_mode,
#' gain_threshold_effective and del_threshold_effective
#' @keywords internal
pdf_cna_plot_settings <- function(report) {
  cna_conf <- report$settings$conf$somatic_cna
  ploidy <- report$settings$conf$sample_properties$tumor_ploidy
  ploidy <- if (is.null(ploidy) || is.na(ploidy) || ploidy == "NA") 2 else
    as.numeric(ploidy)
  amp <- cna_conf$amp_threshold_effective
  if (is.null(amp) || is.na(amp)) amp <- cna_conf$amp_threshold_absolute
  mode <- as.character(cna_conf$threshold_mode)
  relative <- identical(mode, "relative")
  list(
    tumor_ploidy = ploidy,
    amp_threshold_effective = as.numeric(amp),
    threshold_mode = mode,
    gain_threshold_effective = if (relative)
      ploidy * as.numeric(cna_conf$gain_threshold_relative) else
        as.numeric(cna_conf$gain_threshold_absolute),
    del_threshold_effective = if (relative)
      ploidy * as.numeric(cna_conf$del_threshold_relative) else
        as.numeric(cna_conf$del_threshold_absolute))
}
