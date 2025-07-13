#' Plot allele-specific copy number segments
#'
#' Function that plots allele-specific copy number segments
#' (minor + total allele copies)
#'
#' @param chrom_coordinates data frame with assembly-specific chromosome coordinate data (length etc)
#' @param cna_segment data frame with annotated copy number segments
#' @param cna_gene data frame with gene-level copy number data
#'
#' @export
#'
plot_cna_segments <- function(chrom_coordinates = NULL,
                              cna_segment = NULL,
                              cna_gene = NULL){

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


  ## Identify segments that involve oncogene gain or
  ## tumor suppressor loss (homozygous or heterozygous)
  onc_gain_tsg_loss <- cna_gene |>
    dplyr::select(
      c("CHROM", "SEGMENT_START", "SEGMENT_END",
        "VARIANT_CLASS","ONCOGENE", "ONCOGENE_RANK",
        "TUMOR_SUPPRESSOR","TUMOR_SUPPRESSOR_RANK", "SYMBOL")) |>
    dplyr::mutate(CHROM = paste0("chr", .data$CHROM)) |>
    dplyr::filter(
      (.data$ONCOGENE == TRUE &
         .data$VARIANT_CLASS == "gain") |
        (.data$TUMOR_SUPPRESSOR == TRUE &
           (.data$VARIANT_CLASS == "homdel" |
              .data$VARIANT_CLASS == "hetdel")))

  tsg_loss <- data.frame()
  tsg_het_loss <- data.frame()
  onc_gain <- data.frame()

  ## If there are oncogene gains or tumor suppressor losses,
  ## prepare data for plotting
  if(NROW(onc_gain_tsg_loss) > 0){
    onc_gain <- onc_gain_tsg_loss |>
      dplyr::filter(
        .data$ONCOGENE == TRUE &
          .data$VARIANT_CLASS == "gain")

    ## For now, if multiple oncogenes are involved in an amplified segment, we will only
    ## show the top three in the plot (hover)
    if(NROW(onc_gain) > 0){
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

    tsg_loss <- onc_gain_tsg_loss |>
      dplyr::filter(
        .data$TUMOR_SUPPRESSOR == TRUE &
          .data$VARIANT_CLASS == "homdel")

    ## For now, if multiple TSGs are involved in a lost segment, we will only
    ## show the top three in the plot (hover)
    if(NROW(tsg_loss) > 0){
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

    tsg_het_loss <- onc_gain_tsg_loss |>
      dplyr::filter(
        .data$TUMOR_SUPPRESSOR == TRUE &
          .data$VARIANT_CLASS == "hetdel")

    ## For now, if multiple TSGs are involved in a lost segment, we will only
    ## show the top three in the plot (hover)
    if(NROW(tsg_het_loss) > 0){
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
        "CYTOBAND","EVENT_TYPE","CN_MINOR",
        "CN_TOTAL","EVENT_TYPE")) |>
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
    dplyr::mutate(
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


  ## Add information about lost tumor suppressors and gained oncogenes
  ## to SegmentInfo column
  if(NROW(tsg_loss) > 0){
    cna_segments_global <- cna_segments_global |>
      dplyr::left_join(
        tsg_loss,
        by = c("CHROM", "SEGMENT_START", "SEGMENT_END")
      ) |>
      dplyr::mutate(SegmentInfo = dplyr::if_else(
        !is.na(.data$TSG_LOSS),
        paste0(
          .data$SegmentInfo,
          "<br> - Tumor suppressor loss (homozygous del): ",
          .data$TSG_LOSS),
        .data$SegmentInfo))
  }else{
    cna_segments_global$TSG_LOSS <- as.character(NA)
  }

  if(NROW(tsg_het_loss) > 0){
    cna_segments_global <- cna_segments_global |>
      dplyr::left_join(
        tsg_het_loss,
        by = c("CHROM", "SEGMENT_START", "SEGMENT_END")
      ) |>
      dplyr::mutate(SegmentInfo = dplyr::if_else(
        !is.na(.data$TSG_HET_LOSS),
        paste0(.data$SegmentInfo,
               "<br> - Tumor suppressor loss (heterozygous del): ",
               .data$TSG_HET_LOSS),
        .data$SegmentInfo))
  }else{
    cna_segments_global$TSG_HET_LOSS <- as.character(NA)
  }

  if(NROW(onc_gain) > 0){
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
  if(y_max %% 2 == 0){
    y_max_display <- y_max + 2
  }
  if(y_max > 15){
    y_max_display <- ceiling(y_max/5)*5
  }

  y_axis_interval <- 1
  y_breaks <- seq(0, y_max_display, by = y_axis_interval)
  if(y_max > 5 & y_max <= 10){
    y_axis_interval <- 2
    y_breaks <- c(0,1, seq(2, y_max_display, by = y_axis_interval))
  }else{
    if(y_max > 10 & y_max <= 15){
      y_axis_interval <- 2
      y_breaks <- c(
        0,1, seq(2, y_max_display, by = y_axis_interval))
    }else{
      if(y_max > 15 & y_max <= 30){
        y_axis_interval <- 5
        y_breaks <- c(
          0,1,2,4,seq(5, y_max_display, by = y_axis_interval))
      }else{
        if(y_max > 30){
          y_axis_interval <- 10
          y_breaks <- c(
            0,1,2,5, seq(10, y_max_display, by = y_axis_interval))
        }
      }
    }
  }

  ## Make plot - total copy number + minor copy number
  cna_plot <- ggplot2::ggplot(
    cna_segments_global,
    ggplot2::aes(
      x = .data$SegmentStart,
      y = .data$CN_TOTAL, z = .data$SegmentInfo)) +
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
      data = cna_segments_global |>
        dplyr::mutate(Track = "Total copy number"),
      ggplot2::aes(
        x = .data$SegmentStart,
        xend = .data$SegmentEnd,
        y = .data$CN_TOTAL,
        yend = .data$CN_TOTAL,
        colour = .data$Track), linewidth = 1.6) +
    ggplot2::scale_color_manual(
      values = c(`Total copy number` = "black",
                 `Minor copy number` = "firebrick3")) +
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
      breaks = y_breaks) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(
        0, 0, 0, 0))

}


#' Get oncogenic copy number events
#'
#' Function that extracts oncogenic copy number events from
#' a data frame of annotated transcripts (within copy number segments),
#' utilizing oncogene/tumor suppressor status and copy number
#' variant class (gain/loss). This set is used to highlight
#' copy-altered transcripts that are presumably oncogenic, yet
#' not actionable _per se_.
#'
#' @param cna_df_display data frame with transcript annotations per
#' copy number segment
#'
#' @export
#'
get_oncogenic_cna_events <- function(cna_df_display = NULL){

  cna_oncogenic_events <- data.frame()

  assertthat::assert_that(
    !is.null(cna_df_display)
  )
  assertthat::assert_that(
    is.data.frame(cna_df_display))
  if(NROW(cna_df_display) == 0){
    return(cna_oncogenic_events)
  }
  assertable::assert_colnames(
    cna_df_display,
    c("ACTIONABILITY_TIER",
      "ONCOGENE",
      "TUMOR_SUPPRESSOR",
      "VARIANT_CLASS"),
    only_colnames = FALSE,
    quiet = TRUE)

  assertthat::assert_that(
    is.numeric(cna_df_display$ACTIONABILITY_TIER),
    is.logical(cna_df_display$ONCOGENE),
    is.logical(cna_df_display$TUMOR_SUPPRESSOR),
    is.character(cna_df_display$VARIANT_CLASS)
  )

  oncogene_gain_variants <-
    dplyr::filter(
      cna_df_display,
      !is.na(.data$ACTIONABILITY_TIER) &
        .data$ACTIONABILITY_TIER == 3 &
        .data$ONCOGENE == TRUE &
        .data$VARIANT_CLASS == "gain") |>
    dplyr::select(
      dplyr::any_of(
        pcgrr::dt_display$cna_other_oncogenic
      )
    )

  tsgene_loss_variants <-
    dplyr::filter(
      cna_df_display,
      !is.na(.data$ACTIONABILITY_TIER) &
        .data$ACTIONABILITY_TIER == 3 &
        .data$TUMOR_SUPPRESSOR == TRUE &
        .data$VARIANT_CLASS == "homdel") |>
    dplyr::select(
      dplyr::any_of(
        pcgrr::dt_display$cna_other_oncogenic
      )
    )

  tsgene_hetloss_variants <-
    dplyr::filter(
      cna_df_display,
        .data$TUMOR_SUPPRESSOR == TRUE &
        .data$VARIANT_CLASS == "hetdel") |>
    dplyr::select(
      dplyr::any_of(
        pcgrr::dt_display$cna_other_oncogenic
      )
    )

  cna_oncogenic_events <-
    dplyr::bind_rows(
      oncogene_gain_variants,
      tsgene_loss_variants,
      tsgene_hetloss_variants
    ) |>
    dplyr::select(
      dplyr::any_of(
        pcgrr::dt_display$cna_other_oncogenic
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$TISSUE_ASSOC_RANK),
      dplyr::desc(.data$GLOBAL_ASSOC_RANK),
    )

  if("SEGMENT_LENGTH_MB" %in% colnames(cna_oncogenic_events)){
    cna_oncogenic_events <- cna_oncogenic_events |>
      dplyr::mutate(SEGMENT_LENGTH_MB = round(
        .data$SEGMENT_LENGTH_MB, digits = 2))
  }

  return(cna_oncogenic_events)


}
