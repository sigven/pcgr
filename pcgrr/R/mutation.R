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
    font_family = "Helvetica"){

  invisible(assertthat::assert_that(
    !is.null(report),
    msg = "Argument 'report' must be provided"))
  invisible(assertthat::assert_that(
    is.list(report),
    msg = "Argument 'report' must be a valid PCGR report object"))
  invisible(assertthat::assert_that(
    !is.null(report$content$snv_indel$callset$variant),
    msg = "Argument 'report' must contain a valid callset"))

  vaf_dist_tumor <- pcgrr::af_distribution(
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
        x = bin_start,
        y = Count,
        fill = VARIANT_CLASS),
      stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::ylab("Number of variants") +
    ggplot2::xlab("Variant allelic fraction - tumor") +
    ggplot2::scale_fill_manual(
      values = pcgrr::color_palette$tier$values) +
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
        vjust = -2),
      axis.text.y = ggplot2::element_text(
        family = "Helvetica",
        size = font_size),
      axis.title.y = ggplot2::element_text(
        family = "Helvetica",
        size = font_size,
        vjust = 1.5),
      plot.margin = (grid::unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
      legend.text = ggplot2::element_text(
        family =
          font_family,
        size = font_size))

  vaf_plot_plotly <- plotly::ggplotly(vaf_plot)
  vaf_plot_plotly$x$layout$legend$title$text <- ""

  return(vaf_plot_plotly)

}
