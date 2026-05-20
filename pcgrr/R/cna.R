#' Plot allele-specific copy number segments (absolute copies)
#'
#' Function that plots allele-specific copy number segments
#' (minor + total allele copies)
#'
#' @param chrom_coordinates data frame with assembly-specific chromosome coordinate data (length etc)
#' @param cna_segment data frame with annotated copy number segments
#' @param cna_gene data frame with gene-level copy number data
#' @param tumor_ploidy numeric tumor ploidy used as the neutral baseline.
#'   A dotted reference line is always drawn at this value. Default 2.
#' @param amp_threshold_effective numeric effective amplification threshold in
#'   absolute copy number units. A dotted reference line is drawn at this value
#'   when \code{threshold_mode} is "absolute" or "combined". Default 5.
#' @param threshold_mode character thresholding mode applied to all CNA tiers:
#'   "absolute", "relative", or "combined". The amplification threshold line is shown only
#'   when mode is "absolute" or "combined". Default "absolute".
#' @param gain_threshold_effective numeric effective gain threshold in absolute
#'   copy number units. A dotted reference line is drawn at this value. Default 3.
#' @param del_threshold_effective numeric effective heterozygous deletion threshold
#'   in absolute copy number units. A dotted reference line is drawn at this value. Default 1.
#'
#' @export
#'
plot_cna_segments_absolute <- function(
    chrom_coordinates = NULL,
    cna_segment = NULL,
    cna_gene = NULL,
    tumor_ploidy = 2,
    amp_threshold_effective = 5,
    threshold_mode = "absolute",
    gain_threshold_effective = 3,
    del_threshold_effective = 1) {

  ## Validate input
  invisible(assertthat::assert_that(
    !is.null(chrom_coordinates),
    !is.null(cna_segment),
    !is.null(cna_gene)
  ))

  invisible(assertthat::assert_that(
    base::is.data.frame(cna_segment),
    base::is.data.frame(cna_gene),
    base::is.data.frame(chrom_coordinates)
  ))

  ## Check required column names
  assertable::assert_colnames(
    chrom_coordinates,
    c("chrom",
      "genome_start",
      "genome_end",
      "length",
      "centromere_left",
      "centromere_right"),
    only_colnames = F,
    quiet = T
  )

  ## Add centromere and midpoint coordinates to chrom_coordinates
  reference_coordinates <- chrom_coordinates |>
    dplyr::mutate(chrom = paste0("chr", .data$chrom)) |>
    dplyr::mutate(
      centromere_genome_start =
        .data$genome_start + .data$centromere_left,
      centromere_genome_end =
        .data$genome_start + .data$centromere_right,
      midpoint = round(
        (.data$genome_start + .data$genome_end) / 2))


  ## Check required column names of cna_segments data frame
  assertable::assert_colnames(
    cna_segment,
    c("CHROM",
      "SEGMENT_START",
      "SEGMENT_END",
      "CN_MAJOR",
      "CN_MINOR",
      "CN_TOTAL",
      "VARIANT_CLASS",
      "EVENT_TYPE"),
    only_colnames = F,
    quiet = T
  )

  ## Check required column names of cna_genes data frame
  assertable::assert_colnames(
    cna_gene,
    c("CHROM",
      "SEGMENT_START",
      "SEGMENT_END",
      "CN_MAJOR",
      "CN_MINOR",
      "CN_TOTAL",
      "ONCOGENE",
      "ONCOGENE_RANK",
      "TUMOR_SUPPRESSOR",
      "TUMOR_SUPPRESSOR_RANK",
      "SYMBOL",
      "VARIANT_CLASS"),
    only_colnames = F,
    quiet = T
  )


  ## Identify segments that involve oncogene amplification/gain or
  ## tumor suppressor loss (homozygous or heterozygous)
  onc_ampl_tsg_loss <- cna_gene |>
    dplyr::select(
      c("CHROM", "SEGMENT_START", "SEGMENT_END",
        "VARIANT_CLASS","ONCOGENE", "ONCOGENE_RANK",
        "TUMOR_SUPPRESSOR","TUMOR_SUPPRESSOR_RANK", "SYMBOL")) |>
    dplyr::mutate(CHROM = paste0("chr", .data$CHROM)) |>
    dplyr::filter(
      (.data$ONCOGENE == TRUE &
         (.data$VARIANT_CLASS == "amplification" |
            .data$VARIANT_CLASS == "gain")) |
        (.data$TUMOR_SUPPRESSOR == TRUE &
           (.data$VARIANT_CLASS == "homdel" |
              .data$VARIANT_CLASS == "hemdel" |
              .data$VARIANT_CLASS == "hetdel")))

  tsg_loss <- data.frame()
  tsg_het_loss <- data.frame()
  onc_ampl <- data.frame()
  onc_gain <- data.frame()

  ## If there are oncogene gains/amplifications or tumor suppressor losses,
  ## prepare data for plotting
  if (NROW(onc_ampl_tsg_loss) > 0) {
    onc_ampl <- onc_ampl_tsg_loss |>
      dplyr::filter(
        .data$ONCOGENE == TRUE &
          .data$VARIANT_CLASS == "amplification")

    ## For now, if multiple oncogenes are involved in an amplified segment, we will only
    ## show the top three in the plot (hover)
    if (NROW(onc_ampl) > 0) {
      onc_ampl <- onc_ampl |>
        dplyr::arrange(dplyr::desc(.data$ONCOGENE_RANK)) |>
        dplyr::group_by(
          .data$CHROM,
          .data$SEGMENT_START,
          .data$SEGMENT_END) |>
        dplyr::summarise(
          ONC_AMPL = paste(
            utils::head(.data$SYMBOL, 3), collapse = ", "),
          .groups = "drop")
    }

    onc_gain <- onc_ampl_tsg_loss |>
      dplyr::filter(
        .data$ONCOGENE == TRUE &
          .data$VARIANT_CLASS == "gain")

    if (NROW(onc_gain) > 0) {
      onc_gain <- onc_gain |>
        dplyr::arrange(dplyr::desc(.data$ONCOGENE_RANK)) |>
        dplyr::group_by(
          .data$CHROM,
          .data$SEGMENT_START,
          .data$SEGMENT_END) |>
        dplyr::summarise(
          ONC_GAIN = paste(
            utils::head(.data$SYMBOL, 3), collapse = ", "),
          .groups = "drop")
    }

    tsg_loss <- onc_ampl_tsg_loss |>
      dplyr::filter(
        .data$TUMOR_SUPPRESSOR == TRUE &
          (.data$VARIANT_CLASS == "hemdel" |
          .data$VARIANT_CLASS == "homdel"))

    ## For now, if multiple TSGs are involved in a lost segment, we will only
    ## show the top three in the plot (hover)
    if (NROW(tsg_loss) > 0) {
      tsg_loss <- tsg_loss |>
        dplyr::arrange(
          dplyr::desc(.data$TUMOR_SUPPRESSOR_RANK)) |>
        dplyr::group_by(
          .data$CHROM,
          .data$SEGMENT_START,
          .data$SEGMENT_END) |>
        dplyr::summarise(
          TSG_LOSS = paste(
            utils::head(.data$SYMBOL, 3), collapse = ", "),
          .groups = "drop")
    }

    tsg_het_loss <- onc_ampl_tsg_loss |>
      dplyr::filter(
        .data$TUMOR_SUPPRESSOR == TRUE &
          .data$VARIANT_CLASS == "hetdel")

    ## For now, if multiple TSGs are involved in a lost segment, we will only
    ## show the top three in the plot (hover)
    if (NROW(tsg_het_loss) > 0) {
      tsg_het_loss <- tsg_het_loss |>
        dplyr::arrange(
          dplyr::desc(.data$TUMOR_SUPPRESSOR_RANK)) |>
        dplyr::group_by(
          .data$CHROM,
          .data$SEGMENT_START,
          .data$SEGMENT_END) |>
        dplyr::summarise(
          TSG_HET_LOSS = paste(
            utils::head(.data$SYMBOL, 3), collapse = ", "),
          .groups = "drop")
    }

  }

  ## Prepare data for plotting
  ## - pull out core segment elements - ignoring gene/transcript annotations
  ## - add segment size in Mb/Kb
  ## - add segment start and end positions in genome coordinates
  cna_segments_global <- cna_segment |>
    dplyr::select(
      c("CHROM", "SEGMENT_START", "SEGMENT_END",
        "CYTOBAND","CN_MINOR",
        "CN_TOTAL","EVENT_TYPE","VARIANT_CLASS")) |>
    # dplyr::left_join(
    #   dplyr::select(
    #     cna_gene, c("CHROM", "SEGMENT_START", "SEGMENT_END",
    #           "VARIANT_CLASS")),
    #   by = c("CHROM", "SEGMENT_START", "SEGMENT_END")
    # ) |>
    dplyr::mutate(CHROM = paste0("chr", .data$CHROM)) |>
    dplyr::left_join(
      dplyr::select(
        reference_coordinates,
        c("chrom", "genome_start","genome_end")),
      by = c("CHROM" = "chrom"),
    ) |>
    dplyr::mutate(
      segsize = round((.data$SEGMENT_END -
                         .data$SEGMENT_START) / 1000000,
                      digits = 3)) |>
    dplyr::mutate(segsize = dplyr::if_else(
      .data$segsize >= 1,
      paste0(round(.data$segsize, digits = 1), "Mb"),
      paste0(.data$segsize * 1000, "Kb"))) |>

    append_styled_cna_vclass(
      colname = "SegmentClass") |>
    dplyr::mutate(
      # SegmentClass = dplyr::case_when(
      #   .data$VARIANT_CLASS == "amplification"          ~ "Amplification",
      #   .data$VARIANT_CLASS == "gain"                   ~ "Gain",
      #   .data$VARIANT_CLASS == "hetdel"                 ~ "Shallow deletion",
      #   .data$VARIANT_CLASS %in% c("homdel", "hemdel") ~ "Deep deletion",
      #   TRUE                                         ~ "No aberration/neutral"
      # ),
      SegmentStart = .data$genome_start + .data$SEGMENT_START,
      SegmentEnd = .data$genome_start + .data$SEGMENT_END,
      SegmentInfo = paste0(paste(
        .data$CHROM, paste(
          scales::comma(.data$SEGMENT_START),
          scales::comma(.data$SEGMENT_END),
          sep = "-"),
        sep = ":"), " (",.data$segsize,")<br> - ",
        .data$CYTOBAND, " (", .data$EVENT_TYPE,")")) |>
    dplyr::select(
      -c("genome_start","EVENT_TYPE","segsize")) |>
    dplyr::distinct()


  ## Add information about lost tumor suppressors
  ## and amplified oncogenes
  ## to SegmentInfo column
  if (NROW(tsg_loss) > 0) {
    cna_segments_global <- cna_segments_global |>
      dplyr::left_join(
        tsg_loss,
        by = c("CHROM", "SEGMENT_START", "SEGMENT_END")
      ) |>
      dplyr::mutate(SegmentInfo = dplyr::if_else(
        !is.na(.data$TSG_LOSS),
        paste0(
          .data$SegmentInfo,
          "<br> - Tumor suppressor loss (deep del): ",
          .data$TSG_LOSS),
        .data$SegmentInfo))
  }else{
    cna_segments_global$TSG_LOSS <- as.character(NA)
  }

  if (NROW(tsg_het_loss) > 0) {
    cna_segments_global <- cna_segments_global |>
      dplyr::left_join(
        tsg_het_loss,
        by = c("CHROM", "SEGMENT_START", "SEGMENT_END")
      ) |>
      dplyr::mutate(SegmentInfo = dplyr::if_else(
        !is.na(.data$TSG_HET_LOSS),
        paste0(.data$SegmentInfo,
               "<br> - Tumor suppressor loss (shallow del): ",
               .data$TSG_HET_LOSS),
        .data$SegmentInfo))
  }else{
    cna_segments_global$TSG_HET_LOSS <- as.character(NA)
  }

  if (NROW(onc_ampl) > 0) {
    cna_segments_global <- cna_segments_global |>
      dplyr::left_join(
        onc_ampl,
        by = c("CHROM", "SEGMENT_START", "SEGMENT_END")
      ) |>
      dplyr::mutate(SegmentInfo = dplyr::if_else(
        !is.na(.data$ONC_AMPL),
        paste0(.data$SegmentInfo,
               "<br> - Oncogene amplification: ",
               .data$ONC_AMPL),
        .data$SegmentInfo))
  }else{
    cna_segments_global$ONC_AMPL <- as.character(NA)
  }

  if (NROW(onc_gain) > 0) {
    cna_segments_global <- cna_segments_global |>
      dplyr::left_join(
        onc_gain,
        by = c("CHROM", "SEGMENT_START", "SEGMENT_END")
      ) |>
      dplyr::mutate(SegmentInfo = dplyr::if_else(
        !is.na(.data$ONC_GAIN),
        paste0(.data$SegmentInfo,
               "<br> - Oncogene gain: ",
               .data$ONC_GAIN),
        .data$SegmentInfo))
  }else{
    cna_segments_global$ONC_GAIN <- as.character(NA)
  }


  low = min(reference_coordinates$genome_start)
  upp = max(reference_coordinates$genome_end)
  y_max <- max(unique(cna_segments_global$CN_TOTAL))

  y_max_display <- y_max + 1
  if (y_max %% 2 == 0) {
    y_max_display <- y_max + 2
  }
  if (y_max > 15) {
    y_max_display <- ceiling(y_max/5)*5
  }

  ## CNA class colours from palette
  cna_colors_abs <- stats::setNames(
    color_palette$cna_variant_class$values,
    color_palette$cna_variant_class$levels_display)
  amp_color_abs    <- unname(cna_colors_abs["Amplification"])
  gain_color_abs   <- unname(cna_colors_abs["Gain"])
  hetdel_color_abs <- unname(cna_colors_abs["Shallow deletion"])
  homdel_color_abs <- unname(cna_colors_abs["Deep deletion"])

  ## Flag and colour setup for reference lines
  show_amp_line_abs <- threshold_mode %in% c("absolute", "combined")
  neutral_line_label_abs <- paste0(
    "Neutral baseline (Ploidy = ", round(tumor_ploidy, 2), ")")
  amp_line_label_abs <- paste0(
    "Amplification threshold (Total CN \u2265 ", amp_threshold_effective, ")")
  gain_line_label_abs <- paste0(
    "Gain threshold (Total CN \u2265 ", gain_threshold_effective, ")")
  del_line_label_abs <- paste0(
    "Shallow del threshold (Total CN \u2264 ", del_threshold_effective, ")")

  if (show_amp_line_abs) {
    y_max_display <- max(y_max_display, ceiling(amp_threshold_effective) + 1)
  }

  leg_colors_abs <- c(
    "Amplification"                  = amp_color_abs,
    "Gain"                           = gain_color_abs,
    "No aberration/neutral"          = "#868686",
    "Shallow deletion"               = hetdel_color_abs,
    "Deep deletion"                  = homdel_color_abs,
    "Minor copy number"              = "#276419",
    stats::setNames("gray40",         neutral_line_label_abs),
    stats::setNames(gain_color_abs,   gain_line_label_abs),
    stats::setNames(hetdel_color_abs, del_line_label_abs))
  leg_breaks_abs <- c(
    "Amplification",
    "Gain",
    "No aberration/neutral",
    "Shallow deletion",
    "Deep deletion",
    "Minor copy number",
    neutral_line_label_abs,
    gain_line_label_abs,
    del_line_label_abs)
  if (show_amp_line_abs) {
    leg_colors_abs <- c(
      leg_colors_abs,
      stats::setNames(amp_color_abs, amp_line_label_abs))
    leg_breaks_abs <- c(leg_breaks_abs, amp_line_label_abs)
  }

  y_axis_interval <- 1
  y_breaks <- seq(0, y_max_display, by = y_axis_interval)
  if (y_max > 5 & y_max <= 10) {
    y_axis_interval <- 2
    y_breaks <- c(0,1, seq(2, y_max_display, by = y_axis_interval))
  }else{
    if (y_max > 10 & y_max <= 15) {
      y_axis_interval <- 2
      y_breaks <- c(
        0,1, seq(2, y_max_display, by = y_axis_interval))
    }else{
      if (y_max > 15 & y_max <= 30) {
        y_axis_interval <- 5
        y_breaks <- c(
          0,1,2,4,seq(5, y_max_display, by = y_axis_interval))
      }else{
        if (y_max > 30) {
          y_axis_interval <- 10
          y_breaks <- c(
            0,1,2,5, seq(10, y_max_display, by = y_axis_interval))
        }
      }
    }
  }

  ## Build plot — neutral baseline, then gain/del threshold lines
  cna_plot <- ggplot2::ggplot(
    cna_segments_global,
    ggplot2::aes(
      x = .data$SegmentStart,
      y = .data$CN_TOTAL, z = .data$SegmentInfo)) +
    ggplot2::geom_hline(
      data = data.frame(yint = tumor_ploidy,
                        label = neutral_line_label_abs),
      ggplot2::aes(yintercept = .data$yint, colour = .data$label),
      linetype = "dotted", linewidth = 0.6
    ) +
    ggplot2::geom_hline(
      data = data.frame(yint = gain_threshold_effective,
                        label = gain_line_label_abs),
      ggplot2::aes(yintercept = .data$yint, colour = .data$label),
      linetype = "dotted", linewidth = 0.6
    ) +
    ggplot2::geom_hline(
      data = data.frame(yint = del_threshold_effective,
                        label = del_line_label_abs),
      ggplot2::aes(yintercept = .data$yint, colour = .data$label),
      linetype = "dotted", linewidth = 0.6
    )

  if (show_amp_line_abs) {
    cna_plot <- cna_plot +
      ggplot2::geom_hline(
        data = data.frame(yint = amp_threshold_effective,
                          label = amp_line_label_abs),
        ggplot2::aes(yintercept = .data$yint, colour = .data$label),
        linetype = "dotted", linewidth = 0.6
      )
  }

  cna_plot <- cna_plot +
    ggplot2::geom_segment(
      data = cna_segments_global |>
        dplyr::mutate(Track = "Minor copy number"),
      ggplot2::aes(
        x = .data$SegmentStart,
        xend = .data$SegmentEnd,
        y = .data$CN_MINOR,
        yend = .data$CN_MINOR,
        colour = .data$Track), linewidth = 1) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$SegmentStart,
        xend = .data$SegmentEnd,
        y = .data$CN_TOTAL,
        yend = .data$CN_TOTAL,
        colour = .data$SegmentClass), linewidth = 1.6) +
    ggplot2::scale_color_manual(
      breaks = leg_breaks_abs,
      values = leg_colors_abs,
      name   = NULL) +
    ggplot2::geom_vline(
      xintercept = unique(
        c(0,as.vector(cna_segments_global$genome_end),upp)),
      linetype="dotted", colour = "gray") +
    ggplot2::scale_x_continuous(
      breaks = c(0, reference_coordinates$midpoint, upp),
      labels =
        c("", gsub(pattern = 'chr',
                   replacement = '',
                   reference_coordinates$chr), "")) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("Absolute allele counts") +
    ggplot2::scale_y_continuous(
      limits = c(0, max(y_breaks)),
      breaks = y_breaks) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13),
      legend.margin = ggplot2::margin(
        0, 0, 0, 0))

  return(cna_plot)

}


#' Plot copy number segments (relative log2 fold change)
#'
#' Function that plots copy number segments as log2 fold change
#' relative to tumor ploidy, i.e. log2(total_cn / ploidy).
#' Segments with gains (log2FC > 0) and losses (log2FC < 0) are
#' shown in contrasting colors. Two reference lines are drawn:
#' one at y = 0 (neutral/diploid state) and one at
#' log2(amp_threshold_effective / ploidy) marking the amplification cutoff.
#'
#' @param chrom_coordinates data frame with assembly-specific chromosome coordinate data
#' @param cna_segment data frame with annotated copy number segments
#' @param cna_gene data frame with gene-level copy number data
#' @param tumor_ploidy numeric tumor ploidy used as the log2FC baseline (default 2)
#' @param amp_threshold_effective numeric effective amplification threshold
#'   in absolute copy number units (converted to log2FC for the reference line)
#' @param threshold_mode character thresholding mode applied to all CNA tiers:
#'   "absolute", "relative", or "combined". The amplification threshold line is shown only
#'   when mode is "relative" or "combined". Default "absolute".
#' @param gain_threshold_effective numeric effective gain threshold in absolute
#'   copy number units (converted to log2FC for the reference line). Default 3.
#' @param del_threshold_effective numeric effective heterozygous deletion threshold
#'   in absolute copy number units (converted to log2FC for the reference line). Default 1.
#'
#' @export
#'
plot_cna_segments_relative <-
  function(chrom_coordinates = NULL,
           cna_segment = NULL,
           cna_gene = NULL,
           tumor_ploidy = 2,
           amp_threshold_effective = 5,
           threshold_mode = "absolute",
           gain_threshold_effective = 3,
           del_threshold_effective = 1) {

    ## Validate input
    invisible(assertthat::assert_that(
      !is.null(chrom_coordinates),
      !is.null(cna_segment),
      !is.null(cna_gene)
    ))

    invisible(assertthat::assert_that(
      base::is.data.frame(cna_segment),
      base::is.data.frame(cna_gene),
      base::is.data.frame(chrom_coordinates)
    ))

    assertthat::assert_that(
      is.numeric(tumor_ploidy),
      tumor_ploidy > 0,
      is.numeric(amp_threshold_effective),
      amp_threshold_effective > 0
    )

    ## Check required column names
    assertable::assert_colnames(
      chrom_coordinates,
      c("chrom", "genome_start", "genome_end", "length",
        "centromere_left", "centromere_right"),
      only_colnames = F, quiet = T
    )

    assertable::assert_colnames(
      cna_segment,
      c("CHROM", "SEGMENT_START", "SEGMENT_END",
        "CN_TOTAL", "CYTOBAND", "EVENT_TYPE"),
      only_colnames = F, quiet = T
    )

    assertable::assert_colnames(
      cna_gene,
      c("CHROM", "SEGMENT_START", "SEGMENT_END",
        "ONCOGENE", "ONCOGENE_RANK",
        "TUMOR_SUPPRESSOR", "TUMOR_SUPPRESSOR_RANK",
        "SYMBOL", "VARIANT_CLASS"),
      only_colnames = F, quiet = T
    )

    ## Add centromere and midpoint coordinates
    reference_coordinates <- chrom_coordinates |>
      dplyr::mutate(chrom = paste0("chr", .data$chrom)) |>
      dplyr::mutate(
        centromere_genome_start =
          .data$genome_start + .data$centromere_left,
        centromere_genome_end =
          .data$genome_start + .data$centromere_right,
        midpoint = round(
          (.data$genome_start + .data$genome_end) / 2))

    ## Identify oncogene amplifications/gains and TSG losses for hover labels
    onc_ampl_tsg_loss <- cna_gene |>
      dplyr::select(
        c("CHROM", "SEGMENT_START", "SEGMENT_END",
          "VARIANT_CLASS", "ONCOGENE", "ONCOGENE_RANK",
          "TUMOR_SUPPRESSOR", "TUMOR_SUPPRESSOR_RANK", "SYMBOL")) |>
      dplyr::mutate(CHROM = paste0("chr", .data$CHROM)) |>
      dplyr::filter(
        (.data$ONCOGENE == TRUE &
           (.data$VARIANT_CLASS == "amplification" |
              .data$VARIANT_CLASS == "gain")) |
          (.data$TUMOR_SUPPRESSOR == TRUE &
             (.data$VARIANT_CLASS == "homdel" |
                .data$VARIANT_CLASS == "hemdel" |
                .data$VARIANT_CLASS == "hetdel")))

    tsg_loss <- data.frame()
    tsg_het_loss <- data.frame()
    onc_ampl <- data.frame()
    onc_gain <- data.frame()

    if (NROW(onc_ampl_tsg_loss) > 0) {
      onc_ampl <- onc_ampl_tsg_loss |>
        dplyr::filter(
          .data$ONCOGENE == TRUE &
            .data$VARIANT_CLASS == "amplification")
      if (NROW(onc_ampl) > 0) {
        onc_ampl <- onc_ampl |>
          dplyr::arrange(dplyr::desc(.data$ONCOGENE_RANK)) |>
          dplyr::group_by(
            .data$CHROM, .data$SEGMENT_START, .data$SEGMENT_END) |>
          dplyr::summarise(
            ONC_AMPL = paste(utils::head(.data$SYMBOL, 3), collapse = ", "),
            .groups = "drop")
      }

      onc_gain <- onc_ampl_tsg_loss |>
        dplyr::filter(
          .data$ONCOGENE == TRUE &
            .data$VARIANT_CLASS == "gain")
      if (NROW(onc_gain) > 0) {
        onc_gain <- onc_gain |>
          dplyr::arrange(dplyr::desc(.data$ONCOGENE_RANK)) |>
          dplyr::group_by(
            .data$CHROM,
            .data$SEGMENT_START,
            .data$SEGMENT_END) |>
          dplyr::summarise(
            ONC_GAIN = paste(utils::head(.data$SYMBOL, 3), collapse = ", "),
            .groups = "drop")
      }

      tsg_loss <- onc_ampl_tsg_loss |>
        dplyr::filter(
          .data$TUMOR_SUPPRESSOR == TRUE &
            (.data$VARIANT_CLASS == "hemdel" |
               .data$VARIANT_CLASS == "homdel"))
      if (NROW(tsg_loss) > 0) {
        tsg_loss <- tsg_loss |>
          dplyr::arrange(dplyr::desc(.data$TUMOR_SUPPRESSOR_RANK)) |>
          dplyr::group_by(
            .data$CHROM,
            .data$SEGMENT_START,
            .data$SEGMENT_END) |>
          dplyr::summarise(
            TSG_LOSS = paste(utils::head(.data$SYMBOL, 3), collapse = ", "),
            .groups = "drop")
      }

      tsg_het_loss <- onc_ampl_tsg_loss |>
        dplyr::filter(
          .data$TUMOR_SUPPRESSOR == TRUE &
            .data$VARIANT_CLASS == "hetdel")
      if (NROW(tsg_het_loss) > 0) {
        tsg_het_loss <- tsg_het_loss |>
          dplyr::arrange(dplyr::desc(.data$TUMOR_SUPPRESSOR_RANK)) |>
          dplyr::group_by(
            .data$CHROM,
            .data$SEGMENT_START,
            .data$SEGMENT_END) |>
          dplyr::summarise(
            TSG_HET_LOSS = paste(
              utils::head(.data$SYMBOL, 3), collapse = ", "),
            .groups = "drop")
      }
    }

    ## log2FC of the amplification, gain and deletion thresholds for reference lines
    amp_log2fc  <- log2(amp_threshold_effective / tumor_ploidy)
    gain_log2fc <- log2(gain_threshold_effective / tumor_ploidy)
    del_log2fc  <- log2(max(del_threshold_effective, 0.1) / tumor_ploidy)
    show_amp_line_rel <- threshold_mode %in% c("relative", "combined")

    ## Prepare segment data: compute log2(CN_TOTAL / ploidy)
    ## Cap CN_TOTAL at 0.1 to avoid log2(0) = -Inf for homozygous deletions
    cna_segments_global <- cna_segment |>
      dplyr::select(
        c("CHROM",
          "SEGMENT_START",
          "SEGMENT_END",
          "CYTOBAND",
          "EVENT_TYPE",
          "VARIANT_CLASS",
          "CN_TOTAL")) |>
      dplyr::mutate(CHROM = paste0("chr", .data$CHROM)) |>
      dplyr::left_join(
        dplyr::select(
          reference_coordinates,
          c("chrom", "genome_start", "genome_end")),
        by = c("CHROM" = "chrom")
      ) |>
      dplyr::mutate(
        segsize = round(
          (.data$SEGMENT_END - .data$SEGMENT_START) / 1000000,
          digits = 3),
        Log2FC = round(log2(
          pmax(.data$CN_TOTAL, 0.1) / tumor_ploidy), digits = 2)
      ) |>
      append_styled_cna_vclass(
        colname = "SegmentClass") |>
      #   SegmentClass = dplyr::case_when(
      #     .data$VARIANT_CLASS == "amplification"          ~ "Amplification",
      #     .data$VARIANT_CLASS == "gain"                   ~ "Gain",
      #     .data$VARIANT_CLASS == "hetdel"                 ~ "Shallow deletion",
      #     .data$VARIANT_CLASS %in% c("homdel", "hemdel") ~ "Deep deletion",
      #     TRUE                                         ~ "Neutral")
      # ) |>
      dplyr::mutate(segsize = dplyr::if_else(
        .data$segsize >= 1,
        paste0(round(.data$segsize, digits = 1), "Mb"),
        paste0(.data$segsize * 1000, "Kb"))) |>
      dplyr::mutate(
        SegmentStart = .data$genome_start + .data$SEGMENT_START,
        SegmentEnd   = .data$genome_start + .data$SEGMENT_END,
        SegmentInfo  = paste0(
          paste(.data$CHROM,
                paste(scales::comma(.data$SEGMENT_START),
                      scales::comma(.data$SEGMENT_END),
                      sep = "-"),
                sep = ":"),
          " (", .data$segsize, ")<br> - ",
          .data$CYTOBAND, " (", .data$EVENT_TYPE, ")")
        #"<br> - log\u2082FC: ", round(.data$Log2FC, 3))
      ) |>
      dplyr::select(-c("genome_start", "EVENT_TYPE", "segsize")) |>
      dplyr::distinct()

    ## Add oncogene/TSG hover annotations
    if (NROW(tsg_loss) > 0) {
      cna_segments_global <- cna_segments_global |>
        dplyr::left_join(
          tsg_loss,
          by = c("CHROM", "SEGMENT_START", "SEGMENT_END")) |>
        dplyr::mutate(SegmentInfo = dplyr::if_else(
          !is.na(.data$TSG_LOSS),
          paste0(.data$SegmentInfo,
                 "<br> - Tumor suppressor loss (deep del): ",
                 .data$TSG_LOSS),
          .data$SegmentInfo))
    } else {
      cna_segments_global$TSG_LOSS <- as.character(NA)
    }

    if (NROW(tsg_het_loss) > 0) {
      cna_segments_global <- cna_segments_global |>
        dplyr::left_join(
          tsg_het_loss,
          by = c("CHROM", "SEGMENT_START", "SEGMENT_END")) |>
        dplyr::mutate(SegmentInfo = dplyr::if_else(
          !is.na(.data$TSG_HET_LOSS),
          paste0(.data$SegmentInfo,
                 "<br> - Tumor suppressor loss (shallow del): ",
                 .data$TSG_HET_LOSS),
          .data$SegmentInfo))
    } else {
      cna_segments_global$TSG_HET_LOSS <- as.character(NA)
    }

    if (NROW(onc_ampl) > 0) {
      cna_segments_global <- cna_segments_global |>
        dplyr::left_join(
          onc_ampl,
          by = c("CHROM", "SEGMENT_START", "SEGMENT_END")) |>
        dplyr::mutate(SegmentInfo = dplyr::if_else(
          !is.na(.data$ONC_AMPL),
          paste0(.data$SegmentInfo,
                 "<br> - Oncogene amplification: ",
                 .data$ONC_AMPL),
          .data$SegmentInfo))
    } else {
      cna_segments_global$ONC_AMPL <- as.character(NA)
    }

    if (NROW(onc_gain) > 0) {
      cna_segments_global <- cna_segments_global |>
        dplyr::left_join(
          onc_gain,
          by = c("CHROM", "SEGMENT_START", "SEGMENT_END")) |>
        dplyr::mutate(SegmentInfo = dplyr::if_else(
          !is.na(.data$ONC_GAIN),
          paste0(.data$SegmentInfo,
                 "<br> - Oncogene gain: ",
                 .data$ONC_GAIN),
          .data$SegmentInfo))
    } else {
      cna_segments_global$ONC_GAIN <- as.character(NA)
    }

    low <- min(reference_coordinates$genome_start)
    upp <- max(reference_coordinates$genome_end)

    ## Y-axis limits and breaks
    y_min_data <- min(cna_segments_global$Log2FC, na.rm = TRUE)
    y_max_data <- max(cna_segments_global$Log2FC, na.rm = TRUE)
    y_limit_min <- min(-3, floor(y_min_data) - 0.5)
    y_limit_max <- if (show_amp_line_rel) {
      max(amp_log2fc + 0.5, ceiling(y_max_data) + 0.5)
    } else {
      ceiling(y_max_data) + 0.5
    }
    y_range <- y_limit_max - y_limit_min
    y_break_by <- if (y_range > 12) 2 else 1
    y_breaks <- seq(ceiling(y_limit_min), floor(y_limit_max), by = y_break_by)

    ## Legend label strings — ploidy lives here, not on the Y-axis
    neutral_line_label <- paste0(
      "Neutral baseline (Ploidy = ", round(tumor_ploidy, 2), ")")


    cna_colors <- stats::setNames(
      color_palette$cna_variant_class$values,
      color_palette$cna_variant_class$levels
    )
    amp_color    <- unname(cna_colors["amplification"])
    gain_color   <- unname(cna_colors["gain"])
    hetdel_color <- unname(cna_colors["hetdel"])
    homdel_color <- unname(cna_colors["homdel"])

    gain_line_label_rel <- paste0(
      "Gain threshold (Total CN \u2265 ", gain_threshold_effective, ")")
    del_line_label_rel <- paste0(
      "Hetdel threshold (Total CN \u2264 ", del_threshold_effective, ")")

    ## Base colour scale (always shown)
    leg_colors <- c(
      "Amplification"                  = amp_color,
      "Gain"                           = gain_color,
      "Shallow deletion"               = hetdel_color,
      "Deep deletion"                  = homdel_color,
      "No aberration/neutral"          = "#808080",
      stats::setNames("gray40",       neutral_line_label),
      stats::setNames(gain_color,     gain_line_label_rel),
      stats::setNames(hetdel_color,   del_line_label_rel)
    )
    leg_breaks <- c(
      "Amplification",
      "Gain",
      "Shallow deletion",
      "Deep deletion",
      "No aberration/neutral",
      neutral_line_label,
      gain_line_label_rel,
      del_line_label_rel)

    ## Optionally extend with amplification threshold line entry
    if (show_amp_line_rel) {
      amp_line_label <- paste0(
        "Amplification threshold (Total CN \u2265 ",
        amp_threshold_effective, ")")
      leg_colors <- c(
        leg_colors,
        stats::setNames(amp_color, amp_line_label))
      leg_breaks <- c(
        leg_breaks,
        amp_line_label)
    }

    ## Build plot — neutral baseline, gain threshold, and del threshold lines
    cna_plot <- ggplot2::ggplot(
      cna_segments_global,
      ggplot2::aes(
        x = .data$SegmentStart,
        y = .data$Log2FC,
        z = .data$SegmentInfo)
    ) +
      ggplot2::geom_hline(
        data = data.frame(
          yint = 0, label =
            neutral_line_label),
        ggplot2::aes(
          yintercept = .data$yint,
          colour = .data$label),
        linetype = "dotted", linewidth = 0.6
      ) +
      ggplot2::geom_hline(
        data = data.frame(
          yint = gain_log2fc,
          label = gain_line_label_rel),
        ggplot2::aes(
          yintercept = .data$yint,
          colour = .data$label),
        linetype = "dotted", linewidth = 0.6
      ) +
      ggplot2::geom_hline(
        data = data.frame(
          yint = del_log2fc,
          label = del_line_label_rel),
        ggplot2::aes(
          yintercept = .data$yint,
          colour = .data$label),
        linetype = "dotted",
        linewidth = 0.6
      )

    ## Amplification threshold reference line — conditional on mode
    if (show_amp_line_rel) {
      cna_plot <- cna_plot +
        ggplot2::geom_hline(
          data = data.frame(
            yint = amp_log2fc,
            label = amp_line_label),
          ggplot2::aes(
            yintercept = .data$yint,
            colour = .data$label),
          linetype = "dotted",
          linewidth = 0.6
        )
    }

    cna_plot <- cna_plot +
      ## Segments colored by gain/loss direction
      ggplot2::geom_segment(
        ggplot2::aes(
          x      = .data$SegmentStart,
          xend   = .data$SegmentEnd,
          y      = .data$Log2FC,
          yend   = .data$Log2FC,
          colour = .data$SegmentClass),
        linewidth = 1.6
      ) +
      ggplot2::scale_color_manual(
        breaks = leg_breaks,
        values = leg_colors,
        name   = NULL
      ) +
      ## Chromosome boundary lines
      ggplot2::geom_vline(
        xintercept = unique(
          c(0, as.vector(cna_segments_global$genome_end), upp)),
        linetype = "dotted", colour = "gray"
      ) +
      ggplot2::scale_x_continuous(
        breaks = c(0, reference_coordinates$midpoint, upp),
        labels = c(
          "", gsub(pattern = "chr", replacement = "",
                   reference_coordinates$chrom), "")
      ) +
      ggplot2::scale_y_continuous(
        limits = c(y_limit_min, y_limit_max),
        breaks = y_breaks
      ) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Chromosome") +
      ggplot2::ylab("log\u2082(Total CN / ploidy)") +
      ggplot2::theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 13),
        legend.margin = ggplot2::margin(0, 0, 0, 0)
      )

    return(cna_plot)

  }


#' Get oncogenic copy number events
#'
#' Function that extracts oncogenic copy number events from
#' a data frame of annotated transcripts (within copy number segments),
#' utilizing oncogene/tumor suppressor status and copy number
#' variant class (gain/loss). This set is used to highlight
#' any copy-number altered transcripts that are presumably
#' oncogenic (e.g. oncogene amplifications, tumor suppressor deletions)
#'
#' @param cna_df_display data frame with transcript annotations per
#' copy number segment
#'
#' @export
#'
get_oncogenic_cna_events <- function(cna_df_display = NULL) {

  cna_oncogenic_events <- data.frame()

  assertthat::assert_that(
    !is.null(cna_df_display)
  )
  assertthat::assert_that(
    is.data.frame(cna_df_display))
  if (NROW(cna_df_display) == 0) {
    return(cna_oncogenic_events)
  }
  assertable::assert_colnames(
    cna_df_display,
    c("ACTIONABLE_GENE",
      "ONCOGENE",
      "ENTREZGENE",
      "ONCOGENICITY_OKB",
      "TUMOR_SUPPRESSOR",
      "VARIANT_CLASS"),
    only_colnames = FALSE,
    quiet = TRUE)

  assertthat::assert_that(
    is.logical(cna_df_display$ACTIONABLE_GENE),
    is.logical(cna_df_display$ONCOGENE),
    is.character(cna_df_display$ONCOGENICITY_OKB),
    is.logical(cna_df_display$TUMOR_SUPPRESSOR),
    is.character(cna_df_display$VARIANT_CLASS)
  )

  oncogene_ampl_variants <-
    dplyr::filter(
      cna_df_display,
        (.data$ONCOGENE == TRUE |
           .data$ONCOGENICITY_OKB == "Oncogenic" |
           .data$ONCOGENICITY_OKB == "Likely Oncogenic" |
           .data$ACTIONABLE_GENE == TRUE) &
        .data$VARIANT_CLASS == "amplification")

  tsgene_loss_variants <-
    dplyr::filter(
      cna_df_display,
        .data$TUMOR_SUPPRESSOR == TRUE &
        (.data$VARIANT_CLASS == "hemdel" |
          .data$VARIANT_CLASS == "hetdel" |
        .data$VARIANT_CLASS == "homdel"))

  cna_oncogenic_events <-
    dplyr::bind_rows(
      oncogene_ampl_variants,
      tsgene_loss_variants
    ) |>
    dplyr::select(
      dplyr::any_of(
        table_display_cols$cna_other_oncogenic
      )
    ) |>
    append_styled_cna_vclass(
      colname = "VARIANT_CLASS_DISPLAY"
    ) |>
    dplyr::arrange(
      .data$CN_TOTAL,
      dplyr::desc(.data$TISSUE_ASSOC_RANK),
      dplyr::desc(.data$GLOBAL_ASSOC_RANK),
    )

  if ("SEGMENT_LENGTH_MB" %in% colnames(cna_oncogenic_events)) {
    cna_oncogenic_events <- cna_oncogenic_events |>
      dplyr::mutate(SEGMENT_LENGTH_MB = round(
        .data$SEGMENT_LENGTH_MB, digits = 2))
  }

  return(cna_oncogenic_events)


}


#' Build display data for potential two-hit events (nested reactable)
#'
#' For CNA records flagged with somatic and/or germline loss-of-function
#' SNV/InDel candidates, this function parses the comma-separated
#' \code{TWOHIT_CANDIDATE_SOMATIC} / \code{TWOHIT_CANDIDATE_GERMLINE} columns
#' (format: \code{VAR_ID;CONSEQUENCE} per entry) and joins them against the
#' somatic and germline variant callsets to produce two data frames suitable
#' for a nested reactable:
#'
#' \describe{
#'   \item{\code{main}}{One row per two-hit CNA gene: SYMBOL, GENENAME,
#'     VARIANT_CLASS, CN_TOTAL, LOH, and a hidden \code{.row_id} key.}
#'   \item{\code{nested}}{One row per overlapping LoF variant: \code{.row_id}
#'     (FK), ORIGIN (Somatic / Germline), ALTERATION,
#'     CONSEQUENCE, VAF_TUMOR (somatic only), CLASSIFICATION.}
#' }
#'
#' @param cna_variant data frame - gene-level CNA callset
#'   (\code{pcg_report$content$cna$callset$variant})
#' @param snv_somatic data frame - somatic SNV/InDel callset
#'   (\code{pcg_report$content$snv_indel$callset$variant})
#' @param snv_germline data frame - germline classified callset
#'   (\code{pcg_report$content$germline_classified$callset$variant})
#'
#' @return Named list with elements \code{main} and \code{nested} (both
#'   data frames, empty if no two-hit candidates found).
#'
#' @export
#'
build_twohit_display_data <- function(
    cna_variant = NULL,
    snv_somatic = NULL,
    snv_germline = NULL) {

  result <- list(main = data.frame(), nested = data.frame())

  if (is.null(cna_variant) || !is.data.frame(cna_variant) ||
      NROW(cna_variant) == 0) {
    return(result)
  }

  required_cols <- c("TWOHIT_CANDIDATE_SOMATIC",
                     "TWOHIT_CANDIDATE_GERMLINE",
                     "VARIANT_CLASS",
                     "LOH",
                     "SYMBOL",
                     "GENENAME")
  if (!all(required_cols %in% colnames(cna_variant))) {
    return(result)
  }

  ## Filter CNA rows that carry at least one two-hit candidate;
  ## sort so deletion LOH events appear before copy_neutral
  cna_twohit <- cna_variant |>
    append_styled_cna_vclass(
      colname = "VARIANT_CLASS_DISPLAY") |>
    dplyr::filter(
      (!is.na(.data$TWOHIT_CANDIDATE_SOMATIC) &
         .data$TWOHIT_CANDIDATE_SOMATIC != ".") |
      (!is.na(.data$TWOHIT_CANDIDATE_GERMLINE) &
         .data$TWOHIT_CANDIDATE_GERMLINE != ".")
    ) |>
    dplyr::mutate(
      .loh_order = dplyr::case_when(
        .data$LOH == "deletion"     ~ 1L,
        .data$LOH == "copy-neutral" ~ 2L,
        TRUE                        ~ 3L)) |>
    dplyr::arrange(.data$.loh_order) |>
    dplyr::select(-".loh_order") |>
    dplyr::mutate(.row_id = dplyr::row_number())

  if (NROW(cna_twohit) == 0) {
    return(result)
  }

  all_nested <- list()

  ## ---- Somatic two-hit candidates ----
  has_somatic <- !is.null(snv_somatic) &&
    is.data.frame(snv_somatic) &&
    NROW(snv_somatic) > 0 &&
    "VAR_ID" %in% colnames(snv_somatic)

  if (has_somatic) {
    som_cols <- intersect(
      c("VAR_ID", "ALTERATION",
        "VAF_TUMOR", "ONCOGENICITY"),
      colnames(snv_somatic))
    snv_som_sub <- snv_somatic |>
      dplyr::select(dplyr::all_of(som_cols)) |>
      dplyr::mutate(VAR_ID = as.character(.data$VAR_ID))

    somatic_rows <- cna_twohit |>
      dplyr::filter(
        !is.na(.data$TWOHIT_CANDIDATE_SOMATIC) &
          .data$TWOHIT_CANDIDATE_SOMATIC != ".") |>
      dplyr::select(".row_id", "TWOHIT_CANDIDATE_SOMATIC") |>
      tidyr::separate_rows("TWOHIT_CANDIDATE_SOMATIC", sep = ",") |>
      tidyr::separate(
        "TWOHIT_CANDIDATE_SOMATIC",
        into = c("VAR_ID", "CONSEQUENCE"),
        sep = ";",
        fill = "right") |>
      dplyr::mutate(
        VAR_ID = trimws(.data$VAR_ID),
        ORIGIN = "Somatic") |>
      dplyr::left_join(snv_som_sub, by = "VAR_ID")

    if ("VAF_TUMOR" %in% colnames(somatic_rows)) {
      somatic_rows <- somatic_rows |>
        dplyr::mutate(VAF_GENOTYPE = dplyr::if_else(
          !is.na(.data$VAF_TUMOR),
          paste0(round(.data$VAF_TUMOR * 100, 1), "%"),
          NA_character_))
    } else {
      somatic_rows$VAF_GENOTYPE <- NA_character_
    }
    if ("ONCOGENICITY" %in% colnames(somatic_rows)) {
      somatic_rows <- dplyr::rename(somatic_rows, CLASSIFICATION = .data$ONCOGENICITY)
    } else {
      somatic_rows$CLASSIFICATION <- NA_character_
    }
    all_nested <- c(all_nested, list(somatic_rows))
  }

  ## ---- Germline two-hit candidates ----
  has_germline <- !is.null(snv_germline) &&
    is.data.frame(snv_germline) &&
    NROW(snv_germline) > 0 &&
    "VAR_ID" %in% colnames(snv_germline)

  if (has_germline) {
    germ_cols <- intersect(
      c("VAR_ID", "ALTERATION",
        "GENOTYPE", "CLASSIFICATION"),
      colnames(snv_germline))
    snv_germ_sub <- snv_germline |>
      dplyr::select(dplyr::all_of(germ_cols)) |>
      dplyr::mutate(VAR_ID = as.character(.data$VAR_ID))

    germline_rows <- cna_twohit |>
      dplyr::filter(
        !is.na(.data$TWOHIT_CANDIDATE_GERMLINE) &
          .data$TWOHIT_CANDIDATE_GERMLINE != ".") |>
      dplyr::select(".row_id", "TWOHIT_CANDIDATE_GERMLINE") |>
      tidyr::separate_rows("TWOHIT_CANDIDATE_GERMLINE", sep = ",") |>
      tidyr::separate(
        "TWOHIT_CANDIDATE_GERMLINE",
        into = c("VAR_ID", "CONSEQUENCE"),
        sep = ";",
        fill = "right") |>
      dplyr::mutate(
        VAR_ID = trimws(.data$VAR_ID),
        ORIGIN = "Germline") |>
      dplyr::left_join(snv_germ_sub, by = "VAR_ID")

    germline_rows$VAF_GENOTYPE <- if ("GENOTYPE" %in% colnames(germline_rows))
      germline_rows$GENOTYPE else NA_character_
    #if ("FINAL_CLASSIFICATION" %in% colnames(germline_rows)) {
    #  germline_rows <- dplyr::rename(germline_rows,
    #                                 CLASSIFICATION = .data$FINAL_CLASSIFICATION)
    #} else {
    #  germline_rows$CLASSIFICATION <- NA_character_
    #}
    all_nested <- c(all_nested, list(germline_rows))
  }

  if (length(all_nested) == 0) {
    return(result)
  }

  nested_df <- dplyr::bind_rows(all_nested) |>
    dplyr::arrange(.data$.row_id) |>
    dplyr::select(dplyr::any_of(
      c(".row_id", "ORIGIN",
        "ALTERATION",
        "CONSEQUENCE",
        "VAF_GENOTYPE",
        "CLASSIFICATION")))

  ## Build pre-rendered HTML for each main row (used by JS details renderer)
  .origin_badge <- function(origin) {
    if (is.na(origin) || origin == "") return("")
    bg <- if (origin == "Germline") "#6A1B9A" else "#1565C0"
    paste0(
      '<span style="background:', bg, ';color:white;font-size:0.94em;',
      'font-weight:600;padding:2px 7px;border-radius:10px;',
      'white-space:nowrap;">', origin, "</span>")
  }

  .th <- function(label, align = "left", width = "") {
    w <- if (nzchar(width)) paste0("width:", width, ";") else ""
    paste0(
      '<th style="padding:5px 10px;text-align:', align, ";", w,
      'background:#f9f9fb;color:#333;font-weight:600;',
      'font-size:0.94em;text-transform:uppercase;letter-spacing:0.05em;',
      'border-bottom:2px solid #ddd;">',
      label, "</th>")
  }

  .td <- function(val, align = "left", mono = FALSE) {
    #ff <- if (mono) "font-family:monospace;" else ""
    paste0(
      '<td style="padding:5px 10px;text-align:', align, '">',
      if (is.na(val) || val == ".") "&ndash;" else val,
      "</td>")
  }

  header_row <- paste0(
    "<tr>",
    .th("Origin", "center", "90px"),
    .th("Alteration"),
    .th("Consequence"),
    .th("VAF / Genotype", "center", "100px"),
    .th("Classification"),
    "</tr>")

  nested_html_per_row <- lapply(
    split(nested_df, nested_df$.row_id),
    function(rows) {
      data_rows <- apply(rows, 1, function(r) {
        paste0(
          '<tr style="border-bottom:1px solid #eee;">',
          .td(.origin_badge(r[["ORIGIN"]]), align = "center"),
          .td(r[["ALTERATION"]], mono = TRUE),
          .td(r[["CONSEQUENCE"]]),
          .td(r[["VAF_GENOTYPE"]], align = "center"),
          .td(r[["CLASSIFICATION"]]),
          "</tr>")
      })
      paste0(
        '<div style="padding:8px 16px 12px 40px;background:#f9f9fb;">',
        '<table style="border-collapse:collapse;width:100%;font-size:0.94em;">',
        "<thead>", header_row, "</thead>",
        "<tbody>", paste(data_rows, collapse = ""), "</tbody>",
        "</table></div>")
    })

  main_df <- cna_twohit |>
    dplyr::select(dplyr::any_of(
      c(".row_id",
        "SYMBOL",
        "GENENAME",
        "VARIANT_CLASS_DISPLAY",
        "CN_TOTAL",
        "LOH"))) |>
    dplyr::mutate(
      NESTED_HTML = vapply(
        as.character(.data$.row_id),
        function(id) nested_html_per_row[[id]],
        character(1)))

  result$main   <- main_df
  result$nested <- nested_df

  return(result)
}
