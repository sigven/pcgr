
#' Function that assigns evidence items for SNVs/InDels to ACMG tiers 1 and 2
#' @param pcg_report_snv_indel report object for snv/indels
#'
#' @return pcg_report_snv_indel data frame with all report elements
#' @export

assign_tier1_tier2_acmg <- function(pcg_report_snv_indel) {

  ## Assign get evidence items to tier 1 and tier 2
  ## TIER 1: evidence items in specific tumor type and
  ## of high clinical evidence (A_B)
  ## TIER 2: evidence items in other tumor types of high
  ## clinical evidence (A_B) and low clinical evidence in
  ## specific tumor type

  unique_variants_tier1 <- data.frame()
  unique_variants_tier2 <- data.frame()

  ## eitems
  eitems_specific_ttype <-
    pcg_report_snv_indel[["clin_eitem"]][["specific_ttype"]]
  eitems_any_ttype <-
    pcg_report_snv_indel[["clin_eitem"]][["any_ttype"]]
  eitems_other_ttype <-
    pcg_report_snv_indel[["clin_eitem"]][["other_ttype"]]


  for (etype in c("diagnostic", "predictive", "prognostic")) {
    if (nrow(eitems_specific_ttype[[etype]][["A_B"]]) > 0) {
      vars <-
        dplyr::select(eitems_specific_ttype[[etype]][["A_B"]],
                      .data$GENOMIC_CHANGE) |>
        dplyr::distinct()
      unique_variants_tier1 <-
        rbind(unique_variants_tier1, vars) |>
        dplyr::distinct()
    }
  }

  for (etype in c("diagnostic", "predictive", "prognostic")) {
    if (nrow(eitems_any_ttype[[etype]][["A_B"]]) > 0) {
      eitems_other_ttype[[etype]][["A_B"]] <-
        eitems_any_ttype[[etype]][["A_B"]]

      if (nrow(eitems_specific_ttype[[etype]][["A_B"]]) > 0) {

        if (pcgrr::check_common_colnames(
          df1 = eitems_any_ttype[[etype]][["A_B"]],
          df2 = eitems_specific_ttype[[etype]][["A_B"]],
          cnames = c("GENOMIC_CHANGE"))) {

          eitems_other_ttype[[etype]][["A_B"]] <-
            dplyr::anti_join(eitems_any_ttype[[etype]][["A_B"]],
                             eitems_specific_ttype[[etype]][["A_B"]],
                             by = c("GENOMIC_CHANGE"))
        }
      }
      if (nrow(eitems_other_ttype[[etype]][["A_B"]]) > 0) {
        if (nrow(unique_variants_tier1) > 0) {
          if (pcgrr::check_common_colnames(
            df1 = unique_variants_tier1,
            df2 = eitems_other_ttype[[etype]][["A_B"]],
            cnames = c("GENOMIC_CHANGE"))){
            eitems_other_ttype[[etype]][["A_B"]] <-
              dplyr::anti_join(eitems_other_ttype[[etype]][["A_B"]],
                               unique_variants_tier1,
                               by = c("GENOMIC_CHANGE"))
          }
        }
        if (nrow(eitems_other_ttype[[etype]][["A_B"]]) > 0) {
          unique_variants_tier2 <- unique_variants_tier2 |>
            dplyr::bind_rows(
              dplyr::select(eitems_other_ttype[[etype]][["A_B"]],
                            .data$GENOMIC_CHANGE)) |>
            dplyr::distinct()
        }
      }
    }
    if (nrow(eitems_specific_ttype[[etype]][["C_D_E"]]) > 0) {
      if (nrow(unique_variants_tier1) > 0) {
        if (pcgrr::check_common_colnames(
          df1 = unique_variants_tier1,
          df2 = eitems_specific_ttype[[etype]][["C_D_E"]],
          cnames = c("GENOMIC_CHANGE"))) {
          eitems_specific_ttype[[etype]][["C_D_E"]] <-
            dplyr::anti_join(
              eitems_specific_ttype[[etype]][["C_D_E"]],
              unique_variants_tier1, by = c("GENOMIC_CHANGE"))
        }
      }
      if (nrow(eitems_specific_ttype[[etype]][["C_D_E"]]) > 0) {
        unique_variants_tier2 <- unique_variants_tier2 |>
          dplyr::bind_rows(
            dplyr::select(eitems_specific_ttype[[etype]][["C_D_E"]],
                          .data$GENOMIC_CHANGE)) |>
          dplyr::distinct()
      }
    }
  }

  pcg_report_snv_indel[["disp"]][["tier1"]] <-
    unique_variants_tier1
  pcg_report_snv_indel[["disp"]][["tier2"]] <-
    unique_variants_tier2
  pcg_report_snv_indel[["clin_eitem"]][["specific_ttype"]] <-
    eitems_specific_ttype
  pcg_report_snv_indel[["clin_eitem"]][["any_ttype"]] <-
    eitems_any_ttype
  pcg_report_snv_indel[["clin_eitem"]][["other_ttype"]] <-
    eitems_other_ttype

  if (nrow(unique_variants_tier1) > 0) {
    if (pcgrr::check_common_colnames(
      df1 = pcg_report_snv_indel[["variant_set"]][["tier1"]],
      df2 = unique_variants_tier1,
      cnames = c("GENOMIC_CHANGE"))) {
      pcg_report_snv_indel[["variant_set"]][["tier1"]] <-
        dplyr::semi_join(pcg_report_snv_indel[["variant_set"]][["tier1"]],
                         unique_variants_tier1, by = c("GENOMIC_CHANGE"))
    }
    if (pcgrr::check_common_colnames(
      df1 = pcg_report_snv_indel[["variant_set"]][["tier2"]],
      df2 = unique_variants_tier1,
      cnames = c("GENOMIC_CHANGE"))) {
      pcg_report_snv_indel[["variant_set"]][["tier2"]] <-
        dplyr::anti_join(pcg_report_snv_indel[["variant_set"]][["tier2"]],
                         unique_variants_tier1, by = c("GENOMIC_CHANGE"))
    }
  }
  else{
    pcg_report_snv_indel[["variant_set"]][["tier1"]] <- data.frame()
  }
  if (nrow(unique_variants_tier2) == 0) {
    pcg_report_snv_indel[["variant_set"]][["tier2"]] <- data.frame()
  } else{
    if (pcgrr::check_common_colnames(
      df1 = pcg_report_snv_indel[["variant_set"]][["tier2"]],
      df2 = unique_variants_tier2,
      cnames = c("GENOMIC_CHANGE"))) {
      pcg_report_snv_indel[["variant_set"]][["tier2"]] <-
        dplyr::semi_join(pcg_report_snv_indel[["variant_set"]][["tier2"]],
                         unique_variants_tier2, by = c("GENOMIC_CHANGE"))
    }
  }

  return(pcg_report_snv_indel)
}

#' Function that assigns evidence items for SCNAs to ACMG tiers 1 and 2
#' @param pcg_report_cna report object for CNAs
#'
#' @return pcg_report_cna data frame with all report elements
#'
#' @export

assign_tier1_tier2_acmg_cna <- function(pcg_report_cna) {

  ## Assign get evidence items to tier 1 and tier 2
  ## TIER 1: evidence items in specific tumor type and
  ## of high clinical evidence (A_B)
  ## TIER 2: evidence items in other tumor types of high
  ## clinical evidence (A_B) and low clinical evidence in
  ## specific tumor type

  unique_variants_tier1 <- data.frame()
  unique_variants_tier2 <- data.frame()

  ## eitems
  eitems_specific_ttype <- pcg_report_cna[["clin_eitem"]][["specific_ttype"]]
  eitems_any_ttype <- pcg_report_cna[["clin_eitem"]][["any_ttype"]]
  eitems_other_ttype <- pcg_report_cna[["clin_eitem"]][["other_ttype"]]

  for (etype in c("diagnostic", "predictive", "prognostic")) {
    if (nrow(eitems_specific_ttype[[etype]][["A_B"]]) > 0) {

      assertable::assert_colnames(eitems_specific_ttype[[etype]][["A_B"]],
                                  c("SYMBOL", "SEGMENT", "CNA_TYPE"),
                                  only_colnames = F, quiet = T)

      vars <- dplyr::select(eitems_specific_ttype[[etype]][["A_B"]],
                            .data$SYMBOL, .data$SEGMENT, .data$CNA_TYPE) |>
        dplyr::distinct()
      unique_variants_tier1 <- rbind(unique_variants_tier1, vars) |>
        dplyr::distinct()
    }
  }

  for (etype in c("diagnostic", "predictive", "prognostic")) {
    if (nrow(eitems_any_ttype[[etype]][["A_B"]]) > 0) {
      eitems_other_ttype[[etype]][["A_B"]] <-
        eitems_any_ttype[[etype]][["A_B"]]

      if (nrow(eitems_specific_ttype[[etype]][["A_B"]]) > 0) {

        if (pcgrr::check_common_colnames(
          df1 = eitems_any_ttype[[etype]][["A_B"]],
          df2 = eitems_specific_ttype[[etype]][["A_B"]],
          cnames = c("SYMBOL", "SEGMENT", "CNA_TYPE"))) {

          eitems_other_ttype[[etype]][["A_B"]] <-
            dplyr::anti_join(eitems_any_ttype[[etype]][["A_B"]],
                             eitems_specific_ttype[[etype]][["A_B"]],
                             by = c("SYMBOL", "SEGMENT", "CNA_TYPE"))
        }
      }
      if (nrow(eitems_other_ttype[[etype]][["A_B"]]) > 0) {
        if (nrow(unique_variants_tier1) > 0) {
          if (pcgrr::check_common_colnames(
            df1 = unique_variants_tier1,
            df2 = eitems_other_ttype[[etype]][["A_B"]],
            cnames = c("SYMBOL", "SEGMENT", "CNA_TYPE"))) {
            eitems_other_ttype[[etype]][["A_B"]] <-
              dplyr::anti_join(eitems_other_ttype[[etype]][["A_B"]],
                               unique_variants_tier1,
                               by = c("SYMBOL", "SEGMENT", "CNA_TYPE"))
          }
        }
        if (nrow(eitems_other_ttype[[etype]][["A_B"]]) > 0) {

          assertable::assert_colnames(eitems_other_ttype[[etype]][["A_B"]],
                                      c("SYMBOL", "SEGMENT", "CNA_TYPE"),
                                      only_colnames = F, quiet = T)

          unique_variants_tier2 <- unique_variants_tier2 |>
            dplyr::bind_rows(
              dplyr::select(eitems_other_ttype[[etype]][["A_B"]],
                            .data$SYMBOL, .data$SEGMENT, .data$CNA_TYPE)) |>
            dplyr::distinct()
        }
      }
    }
    if (nrow(eitems_specific_ttype[[etype]][["C_D_E"]]) > 0) {
      if (nrow(unique_variants_tier1) > 0) {
        if (pcgrr::check_common_colnames(
          df1 = unique_variants_tier1,
          df2 = eitems_specific_ttype[[etype]][["C_D_E"]],
          cnames = c("SYMBOL", "SEGMENT", "CNA_TYPE"))) {
          eitems_specific_ttype[[etype]][["C_D_E"]] <-
            dplyr::anti_join(
              eitems_specific_ttype[[etype]][["C_D_E"]],
              unique_variants_tier1,
              by = c("SYMBOL", "SEGMENT", "CNA_TYPE"))
        }
      }
      if (nrow(eitems_specific_ttype[[etype]][["C_D_E"]]) > 0) {

        assertable::assert_colnames(eitems_specific_ttype[[etype]][["C_D_E"]],
                                    c("SYMBOL", "SEGMENT", "CNA_TYPE"),
                                    only_colnames = F, quiet = T)

        unique_variants_tier2 <- unique_variants_tier2 |>
          dplyr::bind_rows(
            dplyr::select(eitems_specific_ttype[[etype]][["C_D_E"]],
                          .data$SYMBOL, .data$SEGMENT, .data$CNA_TYPE)) |>
          dplyr::distinct()
      }
    }
  }

  pcg_report_cna[["disp"]][["tier1"]] <- unique_variants_tier1
  pcg_report_cna[["disp"]][["tier2"]] <- unique_variants_tier2
  pcg_report_cna[["clin_eitem"]][["specific_ttype"]] <- eitems_specific_ttype
  pcg_report_cna[["clin_eitem"]][["any_ttype"]] <- eitems_any_ttype
  pcg_report_cna[["clin_eitem"]][["other_ttype"]] <- eitems_other_ttype

  return(pcg_report_cna)

}
