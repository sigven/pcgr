#' Function that predicts MSI status based on fraction of indels among calls
#'
#' @param vcf_data_df data frame with somatic mutations/indels
#' @param pcgr_data object with PCGR datasets
#' @param msi_prediction_model statistical model for MSI prediction
#' @param msi_prediction_dataset underlying dataset from TCGA used for
#' development of statistical classifier
#' @param target_size_mb size of targeted genomic region (coding)
#' @param sample_name name of sample
#' @return msi_data
#'
#' @export
predict_msi_status <- function(vcf_data_df, pcgr_data,
                               msi_prediction_model,
                               msi_prediction_dataset,
                               target_size_mb, sample_name = "Test") {

  mutations_valid <- pcgrr::get_valid_chromosomes(
    vcf_data_df,
    chromosome_column = "CHROM",
    bsg = pcgr_data[["assembly"]][["bsg"]])
  mutations_valid <- mutations_valid %>%
    dplyr::select(.data$CHROM, .data$POS, .data$REF, .data$ALT, .data$CONSEQUENCE, .data$SYMBOL,
                  .data$GENOMIC_CHANGE, .data$VARIANT_CLASS, .data$PROTEIN_DOMAIN,
                  .data$GENE_NAME, .data$PROTEIN_CHANGE, .data$MUTATION_HOTSPOT,
                  .data$CLINVAR, .data$TCGA_FREQUENCY, .data$AF_TUMOR, .data$DP_TUMOR,
                  .data$AF_CONTROL, .data$DP_CONTROL, .data$CALL_CONFIDENCE,
                  .data$SIMPLEREPEATS_HIT, .data$WINMASKER_HIT)

  vcf_df_repeatAnnotated <- mutations_valid %>%
    dplyr::mutate(repeatStatus =
                    dplyr::if_else(.data$SIMPLEREPEATS_HIT == T,
                                   "simpleRepeat", as.character(NA))) %>%
    dplyr::mutate(winMaskStatus =
                    dplyr::if_else(.data$WINMASKER_HIT == T,
                                   "winMaskDust", as.character(NA)))

  msi_stats <- data.frame("sample_name" = sample_name, stringsAsFactors = F)

  msi_stats1 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(.data$repeatStatus) & (.data$VARIANT_CLASS == "insertion" |
                                            .data$VARIANT_CLASS == "deletion")) %>%
    dplyr::summarise(repeat_indels = dplyr::n())

  msi_stats2 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(.data$repeatStatus) & .data$VARIANT_CLASS == "SNV") %>%
    dplyr::summarise(repeat_SNVs = dplyr::n())

  msi_stats3 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(.data$repeatStatus)) %>%
    dplyr::summarise(repeat_indelSNVs = dplyr::n())

  winmask_indels <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(.data$winMaskStatus) &
                    (.data$VARIANT_CLASS == "insertion" |
                       .data$VARIANT_CLASS == "deletion")) %>%
    dplyr::summarise(winmask_indels = dplyr::n())

  winmask_snvs <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(.data$winMaskStatus) & .data$VARIANT_CLASS == "SNV") %>%
    dplyr::summarise(winmask_SNVs = dplyr::n())

  winmask_tot <- vcf_df_repeatAnnotated %>%
    dplyr::filter(!is.na(.data$winMaskStatus)) %>%
    dplyr::summarise(winmask_indelSNVs = dplyr::n())

  msi_stats4 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(is.na(.data$repeatStatus) &
                    (.data$VARIANT_CLASS == "insertion" |
                       .data$VARIANT_CLASS == "deletion")) %>%
    dplyr::summarise(nonRepeat_indels = dplyr::n())

  msi_stats5 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(is.na(.data$repeatStatus) & .data$VARIANT_CLASS == "SNV") %>%
    dplyr::summarise(nonRepeat_SNVs = dplyr::n())

  msi_stats6 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(is.na(.data$repeatStatus)) %>%
    dplyr::summarise(nonRepeat_indelSNVs = dplyr::n())

  msi_stats7 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(.data$VARIANT_CLASS == "insertion" |
                    .data$VARIANT_CLASS == "deletion") %>%
    dplyr::summarise(indels = dplyr::n())

  msi_stats8 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(.data$VARIANT_CLASS == "SNV") %>%
    dplyr::summarise(SNVs = dplyr::n())

  msi_stats9 <- vcf_df_repeatAnnotated %>%
    dplyr::summarise(indelSNVs = dplyr::n())

  msi_stats10 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "MLH1" &
        stringr::str_detect(
          .data$CONSEQUENCE,
          "frameshift_|missense_|splice_|stop_|start_|inframe_")) %>%
    dplyr::summarise(MLH1 = dplyr::n())

  msi_stats11 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "MLH3" &
        stringr::str_detect(
          .data$CONSEQUENCE,
          "frameshift_|missense_|splice_|stop_|start_|inframe_")) %>%
    dplyr::summarise(MLH3 = dplyr::n())

  msi_stats12 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "MSH2" &
        stringr::str_detect(
          .data$CONSEQUENCE,
          "frameshift_|missense_|splice_|stop_|start_|inframe_")) %>%
    dplyr::summarise(MSH2 = dplyr::n())

  msi_stats13 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "MSH3" &
        stringr::str_detect(
          .data$CONSEQUENCE, "frameshift_|missense_|splice_|stop_|start_|frame_")) %>%
    dplyr::summarise(MSH3 = dplyr::n())

  msi_stats14 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "MSH6" &
        stringr::str_detect(
          .data$CONSEQUENCE, "frameshift_|missense_|splice_|stop_|start_|frame_")) %>%
    dplyr::summarise(MSH6 = dplyr::n())

  msi_stats15 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "PMS1" &
        stringr::str_detect(
          .data$CONSEQUENCE, "frameshift_|missense_|splice_|stop_|start_|frame_")) %>%
    dplyr::summarise(PMS1 = dplyr::n())

  msi_stats16 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "PMS2" &
        stringr::str_detect(
          .data$CONSEQUENCE, "frameshift_|missense_|splice_|stop_|start_|frame_")) %>%
    dplyr::summarise(PMS2 = dplyr::n())

  msi_stats17 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "POLE" &
        stringr::str_detect(
          .data$CONSEQUENCE, "frameshift_|missense_|splice_|stop_|start_|frame_")) %>%
    dplyr::summarise(POLE = dplyr::n())

  msi_stats18 <- vcf_df_repeatAnnotated %>%
    dplyr::filter(
      .data$SYMBOL == "POLD1" &
        stringr::str_detect(
          .data$CONSEQUENCE, "frameshift_|missense_|splice_|stop_|start_|frame_")) %>%
    dplyr::summarise(POLD1 = dplyr::n())

  msi_stats1$sample_name <- sample_name
  msi_stats2$sample_name <- sample_name
  msi_stats3$sample_name <- sample_name
  msi_stats4$sample_name <- sample_name
  msi_stats5$sample_name <- sample_name
  msi_stats6$sample_name <- sample_name
  msi_stats7$sample_name <- sample_name
  msi_stats8$sample_name <- sample_name
  msi_stats9$sample_name <- sample_name
  msi_stats10$sample_name <- sample_name
  msi_stats11$sample_name <- sample_name
  msi_stats12$sample_name <- sample_name
  msi_stats13$sample_name <- sample_name
  msi_stats14$sample_name <- sample_name
  msi_stats15$sample_name <- sample_name
  msi_stats16$sample_name <- sample_name
  msi_stats17$sample_name <- sample_name
  msi_stats18$sample_name <- sample_name
  winmask_indels$sample_name <- sample_name
  winmask_snvs$sample_name <- sample_name
  winmask_tot$sample_name <- sample_name

  msi_stats <- msi_stats %>%
    dplyr::left_join(msi_stats1, by = "sample_name") %>%
    dplyr::left_join(msi_stats2, by = "sample_name") %>%
    dplyr::left_join(msi_stats3, by = "sample_name") %>%
    dplyr::left_join(msi_stats4, by = "sample_name") %>%
    dplyr::left_join(msi_stats5, by = "sample_name") %>%
    dplyr::left_join(msi_stats6, by = "sample_name") %>%
    dplyr::left_join(msi_stats7, by = "sample_name") %>%
    dplyr::left_join(msi_stats8, by = "sample_name") %>%
    dplyr::left_join(msi_stats9, by = "sample_name") %>%
    dplyr::left_join(msi_stats10, by = "sample_name") %>%
    dplyr::left_join(msi_stats11, by = "sample_name") %>%
    dplyr::left_join(msi_stats12, by = "sample_name") %>%
    dplyr::left_join(msi_stats13, by = "sample_name") %>%
    dplyr::left_join(msi_stats14, by = "sample_name") %>%
    dplyr::left_join(msi_stats15, by = "sample_name") %>%
    dplyr::left_join(msi_stats16, by = "sample_name") %>%
    dplyr::left_join(msi_stats17, by = "sample_name") %>%
    dplyr::left_join(msi_stats18, by = "sample_name") %>%
    dplyr::left_join(winmask_tot, by = "sample_name") %>%
    dplyr::left_join(winmask_snvs, by = "sample_name") %>%
    dplyr::left_join(winmask_indels, by = "sample_name")

  msi_stats$fracWinMaskIndels <-
    msi_stats$winmask_indels / msi_stats$indels
  msi_stats$fracWinMaskSNVs <-
    msi_stats$winmask_SNVs / msi_stats$SNVs
  msi_stats$fracRepeatIndels <-
    msi_stats$repeat_indels / msi_stats$repeat_indelSNVs
  msi_stats$fracNonRepeatIndels <-
    msi_stats$nonRepeat_indels / msi_stats$nonRepeat_indelSNVs
  msi_stats$fracIndels <-
    msi_stats$indels / msi_stats$indelSNVs
  msi_stats$tmb <- as.numeric(msi_stats$indelSNVs) / target_size_mb
  msi_stats$tmb_indel <- as.numeric(msi_stats$indels) / target_size_mb
  msi_stats$tmb_snv <- as.numeric(msi_stats$SNVs) / target_size_mb
  for (stat in c("fracWinMaskIndels", "fracWinMaskSNVs", "fracRepeatIndels",
                 "fracRepeatIndels", "fracNonRepeatIndels", "fracIndels",
                 "tmb", "tmb_snv", "tmb_indel")) {
    if (nrow(msi_stats[is.na(msi_stats[stat]), ]) > 0) {
      msi_stats[is.na(msi_stats[stat]), ][stat] <- 0
    }
  }

  mmr_pol_df <- mutations_valid %>%
    dplyr::filter(
      stringr::str_detect(.data$SYMBOL,
                          "^(MLH1|MLH3|MSH2|MSH3|MSH6|PMS1|PMS2|POLD1|POLE)$") &
        stringr::str_detect(.data$CONSEQUENCE,
                            "frameshift_|missense_|splice_|stop_|inframe_"))
  mmr_pol_df <- dplyr::select(mmr_pol_df, -c(.data$CHROM, .data$POS, .data$REF, .data$ALT))
  mmr_pol_df <- dplyr::rename(mmr_pol_df, GENE = .data$SYMBOL)
  mmr_pol_df <- mmr_pol_df %>%
    dplyr::select(.data$GENE, .data$CONSEQUENCE, .data$PROTEIN_CHANGE,
                  .data$GENE_NAME, .data$VARIANT_CLASS, .data$PROTEIN_DOMAIN,
                  dplyr::everything())

  msi_predictors <- c("fracWinMaskIndels", "fracWinMaskSNVs",
                      "fracRepeatIndels",
                      "fracNonRepeatIndels",
                      "fracIndels", "MLH1", "MLH3", "MSH2",
                      "MSH3", "MSH6", "PMS1",
                      "PMS2", "POLD1", "POLE", "tmb",
                      "tmb_indel", "tmb_snv")
  msi_class <- stats::predict(msi_prediction_model,
                       dplyr::select(msi_stats, msi_predictors))
  if (msi_class == "MSS") {
    msi_stats$predicted_class <- "MSS (Microsatellite stable)"
    #msi_stats$vb <- "MSI status:\nMSS"
    msi_stats$vb <- "MSS"
  }
  else{
    msi_stats$predicted_class <- "MSI.H (Microsatellite instability - high)"
    #msi_stats$vb <- "MSI status:\nMSI - High"
    msi_stats$vb <- "MSI - High"
  }
  log4r_info(paste0("Predicted MSI status: ",
                           msi_stats$predicted_class))
  log4r_info(paste0("MSI - Indel fraction: ",
                           round(msi_stats$fracNonRepeatIndels, digits = 3)))
  msi_data <- list("mmr_pol_variants" = mmr_pol_df,
                   "msi_stats" = msi_stats,
                   "tcga_dataset" = msi_prediction_dataset)

  return(msi_data)

}

#' Function that generates MSI prediction data for PCGR report
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#'
#' @export
generate_report_data_msi <- function(sample_calls,
                                     pcgr_data,
                                     sample_name,
                                     pcgr_config) {

  pcg_report_msi <- pcgrr::init_report(config =pcgr_config,
                                       class = "msi")
  log4r_info("------")
  log4r_info("Predicting microsatellite instability status")

  msi_sample_calls <- sample_calls %>% dplyr::filter(.data$EXONIC_STATUS == "exonic")
  log4r_info(paste0("n = ",
                           nrow(msi_sample_calls),
                           " exonic variants used for MSI prediction"))
  if (nrow(msi_sample_calls) >= 1) {
    pcg_report_msi[["prediction"]] <-
      pcgrr::predict_msi_status(
        msi_sample_calls, pcgr_data,
        msi_prediction_model = pcgr_data[["msi"]][["model"]][["model"]],
        msi_prediction_dataset = pcgr_data[["msi"]][["model"]][["tcga_dataset"]],
        target_size_mb = pcgr_config$assay_props$target_size_mb,
        sample_name = sample_name)
    pcg_report_msi[["eval"]] <- TRUE
  }
  else{
    log4r_info("Missing variants for MSI prediction")
    pcg_report_msi[["missing_data"]] <- TRUE
  }

  return(pcg_report_msi)
}

#' Function that plots the indel fraction for a given sample and
#' contrasts this with the distribution for MSI-H/MSS samples from TCGA
#'
#' @param tcga_msi_dataset underlying dataset from TCGA used for development
#' of statistical classifier
#' @param indel_fraction fraction of indels of all mutations (SNVs + indels)
#' @return p
#'
#' @export

msi_indel_fraction_plot <- function(tcga_msi_dataset, indel_fraction) {

  color_vec <- utils::head(pcgrr::color_palette[["tier"]][["values"]],2)
  names(color_vec) <- c("MSS", "MSI.H")

  p <- ggplot2::ggplot(data = tcga_msi_dataset) +
    ggplot2::geom_histogram(mapping = ggplot2::aes(x = .data$fracIndels,
                                                   color = .data$MSI_status,
                                                   fill = .data$MSI_status),
                            position = "dodge", binwidth = 0.01) +
    ggplot2::ylab("Number of TCGA samples") +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::scale_color_manual(values = color_vec) +
    ggplot2::theme_classic() +
    ggplot2::xlab("InDel fraction among somatic calls") +
    ggplot2::theme(plot.title =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 16,
                                           hjust = 0.5, face = "bold"),
                   axis.text.x =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 14, face = "bold"),
                   axis.title.x =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 14),
                   legend.title =
                     ggplot2::element_blank(),
                   legend.text =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 14),
                   axis.text.y =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 14),
                   axis.title.y =
                     ggplot2::element_text(family = "Helvetica",
                                           size = 14, vjust = 1.5),
                   plot.margin = (grid::unit(c(0.5, 2, 2, 0.5), "cm"))) +
    ggplot2::geom_vline(xintercept = as.numeric(indel_fraction),
                        size = 1.4, linetype = "dashed")

  return(p)

}
