
#' Function that plots four value boxes with the most
#' important findings in the cancer genome
#'
#' @param pcg_report pcg report with list elements
#' @return p
#'
#' @export

plot_value_boxes <- function(pcg_report) {
  df <- data.frame(
    x = rep(seq(0, 16, 8), 3),
    y = c(rep(1, 3), rep(4.5, 3), rep(8, 3)),
    h = rep(3, 9),
    w = rep(7, 9),
    info = c(pcg_report[["content"]][["value_box"]][["tmb"]],
             pcg_report[["content"]][["value_box"]][["signatures"]],
             pcg_report[["content"]][["value_box"]][["kataegis"]],
             pcg_report[["content"]][["value_box"]][["tier1"]],
             pcg_report[["content"]][["value_box"]][["tier2"]],
             pcg_report[["content"]][["value_box"]][["scna"]],
             pcg_report[["content"]][["value_box"]][["tumor_purity"]],
             pcg_report[["content"]][["value_box"]][["tumor_ploidy"]],
             pcg_report[["content"]][["value_box"]][["msi"]]

             ),
    color = factor(1:9)
  )

  assay_props <-
    pcg_report[["metadata"]][["config"]][["assay_props"]]

  ## color - tumor-control
  color <- rep(pcgrr::color_palette[["tier"]][["values"]][1], 9)
  if (assay_props[["vcf_tumor_only"]] == T) {
    ## color - tumor-only
    color <- rep(pcgrr::color_palette[["report_color"]][["values"]][2], 9)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, height = .data$h, width = .data$w,
                                        label = .data$info, fill = color)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = "white", fontface = "bold", size = 7) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_manual(values = rep(color, 9)) +
    ggplot2::theme_void() +
    ggplot2::guides(fill = F)

  return(p)
}



