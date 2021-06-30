
#' Function that plots four value boxes with the most
#' important findings in the cancer genome
#'
#' @param pcg_report pcg report with list elements
#' @return p
#'
#'

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

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, height = h, width = w,
                                        label = info, fill = color)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = "white", fontface = "bold", size = 7) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_manual(values = rep(color, 9)) +
    ggplot2::theme_void() +
    ggplot2::guides(fill = F)

  return(p)
}


#' Function that generates value box data for PCGR report
#'
#' @param pcg_report object with existing PCGR report data elements
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#'
generate_report_data_value_box <- function(pcg_report,
                                           pcgr_data,
                                           sample_name,
                                           pcgr_config) {

  pcg_report_value_box <- pcgrr::init_report(config = pcgr_config,
                                             class = "value_box")
  pcgrr:::log4r_info("------")
  pcgrr:::log4r_info("Assigning elements to PCGR value boxes")

  if (!pcg_report[["content"]][["snv_indel"]][["eval"]]) {
    return(pcg_report_value_box)
  }

  rep_cont <- pcg_report[["content"]]
  tumor_properties <- pcgr_config$t_props

  sig_contributions <-
    rep_cont[["m_signature_mp"]][["result"]][["contributions"]]

  if (rep_cont[["m_signature_mp"]][["eval"]]) {
    if (!is.null(sig_contributions)) {
      if (nrow(sig_contributions[["per_group"]]) > 0) {
        ranked_groups <- sig_contributions[["per_group"]] %>%
          dplyr::arrange(desc(prop_group))

        dominant_aetiology <- ranked_groups[1, "group"]
        # pcg_report_value_box[["signatures"]] <-
        #   paste0("Dominant signature:\n", dominant_aetiology)
        pcg_report_value_box[["signatures"]] <- dominant_aetiology
      }
    }
  }

  if (rep_cont[['kataegis']][["eval"]]){
    pcg_report_value_box[["kataegis"]] <- "None"
      num_events <- NROW(rep_cont$kataegis$events)
      if(num_events > 0){
        num_events <- NROW(rep_cont$kataegis$events %>%
                             dplyr::filter(confidence == 3))
        # pcg_report_value_box[["kataegis"]] <-
        #   paste0("Kataegis events:\n", num_events)
        pcg_report_value_box[["kataegis"]] <- num_events
      }
  }

  if (rep_cont[["tumor_purity"]][["eval"]]) {
    if (!is.na(tumor_properties[["tumor_purity"]])) {
      # pcg_report_value_box[["tumor_purity"]] <-
      #   paste0("Tumor purity:\n",
      #          rep_cont[["tumor_purity"]][["estimate"]])
      pcg_report_value_box[["tumor_purity"]] <-
        tumor_properties[["tumor_purity"]]
    }
  }

  if (rep_cont[["tumor_ploidy"]][["eval"]]) {
    if (!is.na(tumor_properties[["tumor_ploidy"]])) {
      # pcg_report_value_box[["tumor_ploidy"]] <-
      #   paste0("Tumor ploidy:\n",
      #          rep_cont[["tumor_ploidy"]][["estimate"]])
      pcg_report_value_box[["tumor_ploidy"]] <-
        tumor_properties[["tumor_ploidy"]]
      }
  }

  if (rep_cont[["tmb"]][["eval"]]) {
    if (!is.null(rep_cont[["tmb"]][["v_stat"]])) {
      # pcg_report_value_box[["tmb"]] <-
      #   paste0("TMB:\n",
      #          rep_cont[["tmb"]][["v_stat"]][["tmb_estimate"]],
      #          " mutations/Mb")
      #
      pcg_report_value_box[["tmb"]] <-
        paste0(rep_cont[["tmb"]][["v_stat"]][["tmb_estimate"]],
               " mutations/Mb")
    }
  }
  if (rep_cont[["msi"]][["eval"]]) {
    if (length(rep_cont[["msi"]][["prediction"]]) > 0) {
      pcg_report_value_box[["msi"]] <-
        rep_cont[["msi"]][["prediction"]][["msi_stats"]][["vb"]]
    }
  }
  if (!is.null(rep_cont[["cna"]])) {
    if (rep_cont[["cna"]][["eval"]]) {
      if (nrow(rep_cont[["cna"]][["disp"]][["tier1"]]) > 0) {
        pcg_report_value_box[["scna"]] <-
          # paste0(
          #   "SCNAs:\n",
            paste(
              unique(
                head(
                  rep_cont[["cna"]][["disp"]][["tier1"]]$SYMBOL, 2)
                ),
              collapse = ", ")
            #)
      }
      else{
        pcg_report_value_box[["scna"]] <-
          "None of strong significance"
      }
    }
  }

  if (rep_cont[["snv_indel"]][["eval"]]) {
    if (length(rep_cont[["snv_indel"]][["variant_set"]]) > 0) {
      if (nrow(rep_cont[["snv_indel"]][["variant_set"]][["tier1"]]) > 0) {
        tier1_genes <-
          unique(
            unlist(rep_cont[["snv_indel"]][["variant_set"]][["tier1"]]$SYMBOL)
            )
        pcg_report_value_box[["tier1"]] <-
          paste(head(tier1_genes, 2), collapse = ", ")
        if (length(tier1_genes) > 2) {
          pcg_report_value_box[["tier1"]] <-
            paste(paste(head(tier1_genes, 2), collapse = ", "),
                  paste(tier1_genes[3:min(4, length(tier1_genes))],
                        collapse = ", "),
                  sep = "\n")
          if (length(tier1_genes) > 4) {
            pcg_report_value_box[["tier1"]] <-
              paste0(pcg_report_value_box[["'tier1"]], "++")
          }
        }
        # pcg_report_value_box[["tier1"]] <-
        #   paste0("Tier 1 variants:\n", pcg_report_value_box[["tier1"]])
      }else{
        pcg_report_value_box[["tier1"]] <- "None"
      }
      if (nrow(rep_cont[["snv_indel"]][["variant_set"]][["tier2"]]) > 0) {
        tier2_genes <-
          unique(
            unlist(
              rep_cont[["snv_indel"]][["variant_set"]][["tier2"]]$SYMBOL))
        pcg_report_value_box[["tier2"]] <-
          paste(head(tier2_genes, 2), collapse = ", ")
        if (length(tier2_genes) > 2) {
          pcg_report_value_box[["tier2"]] <-
            paste(paste(head(tier2_genes, 2), collapse = ", "),
                  paste(tier2_genes[3:min(4, length(tier2_genes))],
                        collapse = ", "),
                  sep = "\n")
          if (length(tier2_genes) > 4) {
            pcg_report_value_box[["tier2"]] <-
              paste0(pcg_report_value_box[["tier2"]], "++")
          }
        }
      }else{
        pcg_report_value_box[["tier2"]] <- "None"
      }
    }
  }
  pcg_report_value_box$eval <- TRUE

  return(pcg_report_value_box)
}

