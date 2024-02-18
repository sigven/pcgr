
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
  eitems_query_ttype <-
    pcg_report_snv_indel[["clin_eitem"]][["query_ttype"]]
  eitems_any_ttype <-
    pcg_report_snv_indel[["clin_eitem"]][["any_ttype"]]
  eitems_other_ttype <-
    pcg_report_snv_indel[["clin_eitem"]][["other_ttype"]]


  for (etype in c("diagnostic", "predictive", "prognostic")) {
    if (nrow(eitems_query_ttype[[etype]][["A_B"]]) > 0) {
      vars <-
        dplyr::select(eitems_query_ttype[[etype]][["A_B"]],
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

      if (nrow(eitems_query_ttype[[etype]][["A_B"]]) > 0) {

        if (pcgrr::check_common_colnames(
          df1 = eitems_any_ttype[[etype]][["A_B"]],
          df2 = eitems_query_ttype[[etype]][["A_B"]],
          cnames = c("GENOMIC_CHANGE"))) {

          eitems_other_ttype[[etype]][["A_B"]] <-
            dplyr::anti_join(eitems_any_ttype[[etype]][["A_B"]],
                             eitems_query_ttype[[etype]][["A_B"]],
                             by = c("GENOMIC_CHANGE"))
        }
      }
      if (nrow(eitems_other_ttype[[etype]][["A_B"]]) > 0) {
        if (nrow(unique_variants_tier1) > 0) {
          if (pcgrr::check_common_colnames(
            df1 = unique_variants_tier1,
            df2 = eitems_other_ttype[[etype]][["A_B"]],
            cnames = c("GENOMIC_CHANGE"))) {
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
    if (nrow(eitems_query_ttype[[etype]][["C_D_E"]]) > 0) {
      if (nrow(unique_variants_tier1) > 0) {
        if (pcgrr::check_common_colnames(
          df1 = unique_variants_tier1,
          df2 = eitems_query_ttype[[etype]][["C_D_E"]],
          cnames = c("GENOMIC_CHANGE"))) {
          eitems_query_ttype[[etype]][["C_D_E"]] <-
            dplyr::anti_join(
              eitems_query_ttype[[etype]][["C_D_E"]],
              unique_variants_tier1, by = c("GENOMIC_CHANGE"))
        }
      }
      if (nrow(eitems_query_ttype[[etype]][["C_D_E"]]) > 0) {
        unique_variants_tier2 <- unique_variants_tier2 |>
          dplyr::bind_rows(
            dplyr::select(eitems_query_ttype[[etype]][["C_D_E"]],
                          .data$GENOMIC_CHANGE)) |>
          dplyr::distinct()
      }
    }
  }

  pcg_report_snv_indel[["disp"]][["tier1"]] <-
    unique_variants_tier1
  pcg_report_snv_indel[["disp"]][["tier2"]] <-
    unique_variants_tier2
  pcg_report_snv_indel[["clin_eitem"]][["query_ttype"]] <-
    eitems_query_ttype
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
  eitems_query_ttype <- pcg_report_cna[["clin_eitem"]][["query_ttype"]]
  eitems_any_ttype <- pcg_report_cna[["clin_eitem"]][["any_ttype"]]
  eitems_other_ttype <- pcg_report_cna[["clin_eitem"]][["other_ttype"]]

  for (etype in c("diagnostic", "predictive", "prognostic")) {
    if (nrow(eitems_query_ttype[[etype]][["A_B"]]) > 0) {

      assertable::assert_colnames(eitems_query_ttype[[etype]][["A_B"]],
                                  c("SYMBOL", "SEGMENT", "CNA_TYPE"),
                                  only_colnames = F, quiet = T)

      vars <- dplyr::select(eitems_query_ttype[[etype]][["A_B"]],
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

      if (nrow(eitems_query_ttype[[etype]][["A_B"]]) > 0) {

        if (pcgrr::check_common_colnames(
          df1 = eitems_any_ttype[[etype]][["A_B"]],
          df2 = eitems_query_ttype[[etype]][["A_B"]],
          cnames = c("SYMBOL", "SEGMENT", "CNA_TYPE"))) {

          eitems_other_ttype[[etype]][["A_B"]] <-
            dplyr::anti_join(eitems_any_ttype[[etype]][["A_B"]],
                             eitems_query_ttype[[etype]][["A_B"]],
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
    if (nrow(eitems_query_ttype[[etype]][["C_D_E"]]) > 0) {
      if (nrow(unique_variants_tier1) > 0) {
        if (pcgrr::check_common_colnames(
          df1 = unique_variants_tier1,
          df2 = eitems_query_ttype[[etype]][["C_D_E"]],
          cnames = c("SYMBOL", "SEGMENT", "CNA_TYPE"))) {
          eitems_query_ttype[[etype]][["C_D_E"]] <-
            dplyr::anti_join(
              eitems_query_ttype[[etype]][["C_D_E"]],
              unique_variants_tier1,
              by = c("SYMBOL", "SEGMENT", "CNA_TYPE"))
        }
      }
      if (nrow(eitems_query_ttype[[etype]][["C_D_E"]]) > 0) {

        assertable::assert_colnames(eitems_query_ttype[[etype]][["C_D_E"]],
                                    c("SYMBOL", "SEGMENT", "CNA_TYPE"),
                                    only_colnames = F, quiet = T)

        unique_variants_tier2 <- unique_variants_tier2 |>
          dplyr::bind_rows(
            dplyr::select(eitems_query_ttype[[etype]][["C_D_E"]],
                          .data$SYMBOL, .data$SEGMENT, .data$CNA_TYPE)) |>
          dplyr::distinct()
      }
    }
  }

  pcg_report_cna[["disp"]][["tier1"]] <- unique_variants_tier1
  pcg_report_cna[["disp"]][["tier2"]] <- unique_variants_tier2
  pcg_report_cna[["clin_eitem"]][["query_ttype"]] <- eitems_query_ttype
  pcg_report_cna[["clin_eitem"]][["any_ttype"]] <- eitems_any_ttype
  pcg_report_cna[["clin_eitem"]][["other_ttype"]] <- eitems_other_ttype

  return(pcg_report_cna)

}

#' Function that assigns tier classifications to somatic CNA segments and
#' SNVs/InDels, based on the presence of biomarker evidence found in
#' the variant set
#'
#' @param vartype variant type ('snv_indel' or 'cna')
#' @param primary_site primary tumor site
#' @param variants_df data frame with variants (SNVs/InDels or CNAs)
#' @param biomarker_items data frame with biomarker evidence items
#'
#' @export
assign_acmg_tiers <- function(
    vartype = "snv_indel",
    primary_site = "Any",
    variants_df = NULL,
    biomarker_items = NULL) {

  invisible(assertthat::assert_that(
    is.data.frame(variants_df),
    msg = paste0("Argument variants_df needs be of type data.frame")))
  assertable::assert_colnames(
    variants_df, c("TUMOR_SUPPRESSOR",
                   "VAR_ID",
                   "VARIANT_CLASS",
                   "ONCOGENE",
                   "ENTREZGENE"),
    only_colnames = F, quiet = T)
  invisible(assertthat::assert_that(
    is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))
  assertable::assert_colnames(
    biomarker_items,
    c("VAR_ID",
      "ENTREZGENE",
      "BM_EVIDENCE_LEVEL",
      "BM_PRIMARY_SITE"),
    only_colnames = F, quiet = T)

  results_acmg <- list()
  tier_classification <- data.frame()

  if (NROW(biomarker_items) > 0) {
    tier_classification <-
      biomarker_items |>
      #results[['biomarker_evidence']][['items']] |>
      dplyr::select(
        c("VAR_ID",
          "VARIANT_CLASS",
          "ENTREZGENE",
          "BM_EVIDENCE_LEVEL",
          "BM_PRIMARY_SITE")) |>
      dplyr::distinct() |>
      dplyr::mutate(ACMG_AMP_TIER = dplyr::case_when(
        .data$BM_PRIMARY_SITE == primary_site &
          primary_site != "Any" &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, "^(A|B)"
          ) ~ as.integer(1),
        .data$BM_PRIMARY_SITE != primary_site &
          #primary_site != "Any" &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, "^(A|B)"
          ) ~ as.integer(2),
        .data$BM_PRIMARY_SITE == primary_site &
          primary_site != "Any" &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, "^(C|D|E)"
          ) ~ as.integer(2),
        TRUE ~ as.integer(100)
      )) |>
      dplyr::group_by(
        .data$VAR_ID,
        .data$ENTREZGENE,
        .data$VARIANT_CLASS) |>
        #c("VAR_ID", "ENTREZGENE", "VARIANT_CLASS")) |>
      dplyr::summarise(
        ACMG_AMP_TIER = min(.data$ACMG_AMP_TIER, na.rm = T),
        .groups = "drop") |>
      dplyr::mutate(ACMG_AMP_TIER = dplyr::if_else(
        .data$ACMG_AMP_TIER == 100,
        as.integer(NA),
        as.integer(.data$ACMG_AMP_TIER)
      ))

    if (vartype == 'snv_indel' &
       "CODING_STATUS" %in% colnames(variants_df)) {

      variants_df <- variants_df |>
        dplyr::left_join(
          tier_classification,
          by = c("VAR_ID","ENTREZGENE","VARIANT_CLASS")) |>

        dplyr::mutate(ACMG_TIER2 = dplyr::if_else(
          (!is.na(.data$TUMOR_SUPPRESSOR) &
             .data$TUMOR_SUPPRESSOR == TRUE) |
            (!is.na(.data$ONCOGENE) &
               .data$ONCOGENE == TRUE) &
            .data$CODING_STATUS == "coding",
          as.integer(3),
          as.integer(NA)
        )) |>
        dplyr::mutate(ACMG_TIER2 = dplyr::if_else(
          is.na(.data$ACMG_TIER2) |
            (!is.na(.data$ACMG_TIER2) &
               .data$ACMG_TIER2 != 3) &
            .data$CODING_STATUS == "coding",
          as.integer(4),
          as.integer(.data$ACMG_TIER2)
        )) |>
        dplyr::mutate(ACMG_AMP_TIER = dplyr::if_else(
          .data$CODING_STATUS == "noncoding",
          as.integer(5),
          as.integer(.data$ACMG_AMP_TIER)
        )) |>
        dplyr::mutate(ACMG_AMP_TIER = dplyr::case_when(
          is.na(.data$ACMG_AMP_TIER) &
            !is.na(.data$ACMG_TIER2) ~ .data$ACMG_TIER2,
          TRUE ~ as.integer(.data$ACMG_AMP_TIER)
        )) |>
        dplyr::select(-c("ACMG_TIER2")) |>
        dplyr::arrange(.data$ACMG_AMP_TIER)
    }else{

      if (vartype == 'cna') {

        variants_df <- variants_df |>
          dplyr::left_join(
            tier_classification,
            by = c("VAR_ID",
                   "ENTREZGENE",
                   "VARIANT_CLASS")) |>
          dplyr::mutate(ACMG_TIER2 = dplyr::if_else(
            (!is.na(.data$TUMOR_SUPPRESSOR) &
               .data$TUMOR_SUPPRESSOR == TRUE &
               .data$VARIANT_CLASS == "homdel") |
              (!is.na(.data$ONCOGENE) &
                 .data$ONCOGENE == TRUE &
                 .data$VARIANT_CLASS == "gain"),
            as.integer(3),
            as.integer(.data$ACMG_AMP_TIER)
          )) |>
          dplyr::mutate(ACMG_AMP_TIER = dplyr::case_when(
            is.na(.data$ACMG_AMP_TIER) &
              !is.na(.data$ACMG_TIER2) ~ .data$ACMG_TIER2,
            TRUE ~ as.integer(.data$ACMG_AMP_TIER)
          )) |>
          dplyr::select(-c("ACMG_TIER2")) |>
          dplyr::arrange(.data$ACMG_AMP_TIER) |>
          dplyr::distinct()

      }
    }
  }
  else{
    if (vartype == 'snv_indel') {
      variants_df <- variants_df |>
        dplyr::mutate(ACMG_AMP_TIER = dplyr::if_else(
          (!is.na(.data$TUMOR_SUPPRESSOR) &
             .data$TUMOR_SUPPRESSOR == TRUE) |
            (!is.na(.data$ONCOGENE) &
               .data$ONCOGENE == TRUE) &
            .data$CODING_STATUS == "coding",
          as.integer(3),
          as.integer(NA)
        )) |>
        dplyr::mutate(ACMG_AMP_TIER = dplyr::if_else(
          is.na(.data$ACMG_AMP_TIER) &
            .data$CODING_STATUS == "coding",
          as.integer(4),
          as.integer(.data$ACMG_AMP_TIER)
        )) |>
        dplyr::mutate(ACMG_AMP_TIER = dplyr::if_else(
          .data$CODING_STATUS == "noncoding",
          as.integer(5),
          as.integer(.data$ACMG_AMP_TIER)
        )) |>
        dplyr::arrange(.data$ACMG_AMP_TIER) |>
        dplyr::distinct()
    }
    if (vartype == 'cna') {

      variants_df <- variants_df |>
        dplyr::mutate(ACMG_AMP_TIER = dplyr::if_else(
          (!is.na(.data$TUMOR_SUPPRESSOR) &
             .data$TUMOR_SUPPRESSOR == TRUE &
             .data$VARIANT_CLASS == "homdel") |
            (!is.na(.data$ONCOGENE) &
               .data$ONCOGENE == TRUE &
               .data$VARIANT_CLASS == "gain"),
          as.integer(3),
          as.integer(NA)
        )) |>
        dplyr::distinct() |>
        dplyr::arrange(.data$ACMG_AMP_TIER)
    }
  }

  biomarker_items <- biomarker_items |>
    dplyr::left_join(
      dplyr::select(
        variants_df,
        c("VAR_ID",
        "ENTREZGENE",
        "ACMG_AMP_TIER")
      ),
      by = c("VAR_ID","ENTREZGENE")
    ) |>
    dplyr::mutate(ACMG_AMP_TIER = dplyr::case_when(
      .data$BM_PRIMARY_SITE == primary_site &
        primary_site != "Any" &
      as.integer(.data$ACMG_AMP_TIER) == 1 &
        stringr::str_detect(
          .data$BM_EVIDENCE_LEVEL,"^(C|D|E)"
        ) ~ as.integer(NA),
      .data$BM_PRIMARY_SITE != primary_site &
        primary_site != "Any" &
        .data$ACMG_AMP_TIER == 2 &
        stringr::str_detect(
          .data$BM_EVIDENCE_LEVEL,"^(C|D|E)"
        ) ~ as.integer(NA),
      .data$BM_PRIMARY_SITE != primary_site &
        primary_site != "Any" &
        .data$ACMG_AMP_TIER == 1 ~ as.integer(NA),
      TRUE ~ as.integer(.data$ACMG_AMP_TIER)
    )) |>
    dplyr::arrange(.data$ACMG_AMP_TIER,
                   .data$BM_EVIDENCE_LEVEL,
                   dplyr::desc(.data$BM_RATING)) |>
    dplyr::distinct()

  results_acmg[['variant']] <- variants_df |>
    dplyr::rename(TIER = .data$ACMG_AMP_TIER) |>
    dplyr::mutate(TIER_GUIDELINE = "ACMG_AMP")

  results_acmg[['biomarker_evidence']][['items']] <-
    biomarker_items |>
    dplyr::rename(TIER = .data$ACMG_AMP_TIER) |>
    dplyr::mutate(TIER_GUIDELINE = "ACMG_AMP")

  results_acmg[['biomarker_evidence']][['tier_classification']] <-
    tier_classification |>
    dplyr::rename(TIER = .data$ACMG_AMP_TIER) |>
    dplyr::mutate(TIER_GUIDELINE = "ACMG_AMP")

  return(results_acmg)

}
