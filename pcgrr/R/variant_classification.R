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
assign_amp_asco_tiers <- function(
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

  results_amp_asco <- list()
  tier_classification <- data.frame()
  biomarker_items_hires <- data.frame()

  if (NROW(biomarker_items) > 0) {

    assertable::assert_colnames(
      biomarker_items,
      c("VAR_ID",
        "ENTREZGENE",
        "BM_EVIDENCE_LEVEL",
        "BM_PRIMARY_SITE"),
      only_colnames = F, quiet = T)

    ## do not consider biomarkers that match variant properties
    ## associated with non-principal (non-canonical) transcripts
    biomarker_items_hires <- biomarker_items |>
      dplyr::filter(
        vartype == "cna" |
          (vartype == "snv_indel" &
             (.data$BM_RESOLUTION != "exon_nonprincipal" &
                .data$BM_RESOLUTION != "hgvsp_nonprincipal")
          )
      )
  }

  if (NROW(biomarker_items_hires) > 0) {
    tier_classification <-
      biomarker_items_hires |>
      dplyr::select(
        c("VAR_ID",
          "VARIANT_CLASS",
          "ENTREZGENE",
          "BM_EVIDENCE_LEVEL",
          "BM_PRIMARY_SITE")) |>
      dplyr::distinct() |>
      dplyr::mutate(AMP_ASCO_TIER = dplyr::case_when(

        # Biomarker site matches primary site of query tumor - strong evidence
        # TIER 1
        .data$BM_PRIMARY_SITE == primary_site &
          primary_site != "Any" &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, "^(A|B)"
          ) ~ as.integer(1),

        ## Biomarker site does not match primary site of query tumor - strong evidence, OR
        ## If tumor site is not given (primary_site = 'Any'), consider any evidence
        ## as tier 2 (weak and strong)

        # TIER 2
        (.data$BM_PRIMARY_SITE != primary_site &
          #primary_site != "Any" &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, "^(A|B)"
          )) |
          primary_site == "Any" ~ as.integer(2),

        ## Biomarker site matches primary site of query tumor - weak evidence

        # TIER 2
        .data$BM_PRIMARY_SITE == primary_site &
          primary_site != "Any" &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, "^(C|D|E)"
          ) ~ as.integer(2),
        TRUE ~ as.integer(100)
      )) |>

      ## Get top tier level for any given variant
      dplyr::group_by(
        .data$VAR_ID,
        .data$ENTREZGENE,
        .data$VARIANT_CLASS) |>
        #c("VAR_ID", "ENTREZGENE", "VARIANT_CLASS")) |>
      dplyr::summarise(
        AMP_ASCO_TIER = min(
          .data$AMP_ASCO_TIER, na.rm = T),
        .groups = "drop") |>
      dplyr::mutate(AMP_ASCO_TIER = dplyr::if_else(
        .data$AMP_ASCO_TIER == 100,
        as.integer(NA),
        as.integer(.data$AMP_ASCO_TIER)
      ))

    if (vartype == 'snv_indel' &
       "CODING_STATUS" %in% colnames(variants_df)) {

      variants_df <- variants_df |>
        dplyr::left_join(
          tier_classification,
          by = c("VAR_ID","ENTREZGENE","VARIANT_CLASS")) |>

        ## Tumor suppressor/oncogene mutations (tier 3)
        ## - low MAF (< 0.001)
        dplyr::mutate(AMP_ASCO_TIER_OTHER_VARS = dplyr::if_else(
          ((!is.na(.data$TUMOR_SUPPRESSOR) &
             .data$TUMOR_SUPPRESSOR == TRUE) |
            (!is.na(.data$ONCOGENE) &
               .data$ONCOGENE == TRUE)) &
            (is.na(.data$gnomADe_AF) |
               .data$gnomADe_AF < 0.001) &
            .data$CODING_STATUS == "coding",
          as.integer(3),
          as.integer(NA)
        )) |>

        ## Tier 4 - as specified in AMP/ASCO guidelines are
        ## not considered here

        ## Tier 5 - other coding mutations
        dplyr::mutate(AMP_ASCO_TIER_OTHER_VARS = dplyr::if_else(
          is.na(.data$AMP_ASCO_TIER_OTHER_VARS) &
            .data$CODING_STATUS == "coding",
          as.integer(5),
          as.integer(.data$AMP_ASCO_TIER_OTHER_VARS)
        )) |>

        ## Tier 6 - other non-coding mutations
        dplyr::mutate(AMP_ASCO_TIER = dplyr::if_else(
          is.na(.data$AMP_ASCO_TIER) &
            .data$CODING_STATUS == "noncoding",
          as.integer(6),
          as.integer(.data$AMP_ASCO_TIER)
        )) |>
        dplyr::mutate(AMP_ASCO_TIER = dplyr::case_when(
          is.na(.data$AMP_ASCO_TIER) &
            !is.na(.data$AMP_ASCO_TIER_OTHER_VARS) ~ .data$AMP_ASCO_TIER_OTHER_VARS,
          TRUE ~ as.integer(.data$AMP_ASCO_TIER)
        )) |>
        dplyr::select(-c("AMP_ASCO_TIER_OTHER_VARS")) |>
        dplyr::arrange(.data$AMP_ASCO_TIER)
    }else{

      if (vartype == 'cna') {

        variants_df <- variants_df |>
          dplyr::left_join(
            tier_classification,
            by = c("VAR_ID",
                   "ENTREZGENE",
                   "VARIANT_CLASS")) |>
          dplyr::mutate(AMP_ASCO_TIER_OTHER_VARS = dplyr::if_else(
            (!is.na(.data$TUMOR_SUPPRESSOR) &
               .data$TUMOR_SUPPRESSOR == TRUE &
               (.data$VARIANT_CLASS == "homdel" |
                  .data$VARIANT_CLASS == "hetdel")) |
              (!is.na(.data$ONCOGENE) &
                 .data$ONCOGENE == TRUE &
                 .data$VARIANT_CLASS == "gain"),
            as.integer(3),
            as.integer(.data$AMP_ASCO_TIER)
          )) |>
          dplyr::mutate(AMP_ASCO_TIER = dplyr::case_when(
            is.na(.data$AMP_ASCO_TIER) &
              !is.na(.data$AMP_ASCO_TIER_OTHER_VARS) ~ .data$AMP_ASCO_TIER_OTHER_VARS,
            TRUE ~ as.integer(.data$AMP_ASCO_TIER)
          )) |>
          dplyr::select(-c("AMP_ASCO_TIER_OTHER_VARS")) |>
          dplyr::arrange(.data$AMP_ASCO_TIER) |>
          dplyr::distinct()

      }
    }
  }else{
    if (vartype == 'snv_indel') {
      variants_df <- variants_df |>
        dplyr::mutate(AMP_ASCO_TIER = dplyr::if_else(
          ((!is.na(.data$TUMOR_SUPPRESSOR) &
             .data$TUMOR_SUPPRESSOR == TRUE) |
            (!is.na(.data$ONCOGENE) &
               .data$ONCOGENE == TRUE)) &
            (is.na(.data$gnomADe_AF) |
            .data$gnomADe_AF < 0.001) &
            .data$CODING_STATUS == "coding",
          as.integer(3),
          as.integer(NA)
        )) |>
        dplyr::mutate(AMP_ASCO_TIER = dplyr::if_else(
          is.na(.data$AMP_ASCO_TIER) &
            .data$CODING_STATUS == "coding",
          as.integer(5),
          as.integer(.data$AMP_ASCO_TIER)
        )) |>
        dplyr::mutate(AMP_ASCO_TIER = dplyr::if_else(
          is.na(.data$AMP_ASCO_TIER) &
            .data$CODING_STATUS == "noncoding",
          as.integer(6),
          as.integer(.data$AMP_ASCO_TIER)
        )) |>
        dplyr::arrange(.data$AMP_ASCO_TIER) |>
        dplyr::distinct()
    }else{
      if (vartype == 'cna') {

        variants_df <- variants_df |>
          dplyr::mutate(AMP_ASCO_TIER = dplyr::if_else(
            (!is.na(.data$TUMOR_SUPPRESSOR) &
               .data$TUMOR_SUPPRESSOR == TRUE &
               (.data$VARIANT_CLASS == "homdel" |
                  .data$VARIANT_CLASS == "hetdel")) |
              (!is.na(.data$ONCOGENE) &
                 .data$ONCOGENE == TRUE &
                 .data$VARIANT_CLASS == "gain"),
            as.integer(3),
            as.integer(NA)
          )) |>
          dplyr::distinct() |>
          dplyr::arrange(.data$AMP_ASCO_TIER)
      }
    }
  }

  if(NROW(biomarker_items_hires) > 0){
    biomarker_items_hires <- biomarker_items_hires |>
      dplyr::left_join(
        dplyr::select(
          variants_df,
          c("VAR_ID",
          "ENTREZGENE",
          "AMP_ASCO_TIER")
        ),
        by = c("VAR_ID",
               "ENTREZGENE")
      ) |>
      dplyr::mutate(AMP_ASCO_TIER = dplyr::case_when(

        ## update variant tier status for individual evidence items
        .data$BM_PRIMARY_SITE == primary_site &
          primary_site != "Any" &
        as.integer(.data$AMP_ASCO_TIER) == 1 &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL,"^(C|D|E)"
          ) ~ as.integer(NA),
        .data$BM_PRIMARY_SITE != primary_site &
          primary_site != "Any" &
          .data$AMP_ASCO_TIER == 2 &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL,"^(C|D|E)"
          ) ~ as.integer(NA),
        .data$BM_PRIMARY_SITE != primary_site &
          primary_site != "Any" &
          .data$AMP_ASCO_TIER == 1 ~ as.integer(NA),
        TRUE ~ as.integer(.data$AMP_ASCO_TIER)
      )) |>
      dplyr::arrange(
        .data$AMP_ASCO_TIER,
        .data$BM_EVIDENCE_LEVEL,
        dplyr::desc(.data$BM_RATING)) |>
      dplyr::distinct()
  }

  results_amp_asco[['biomarker_evidence']][['tier_classification']] <- data.frame()
  results_amp_asco[['biomarker_evidence']][['items']] <- data.frame()
  results_amp_asco[['biomarker_evidence']][['items_all']] <- data.frame()
  results_amp_asco[['variant']] <- data.frame()

  results_amp_asco[['variant']] <- variants_df |>
    dplyr::rename(ACTIONABILITY_TIER = .data$AMP_ASCO_TIER) |>
    dplyr::mutate(ACTIONABILITY_FRAMEWORK = "AMP_ASCO_CAP") |>
    dplyr::mutate(ACTIONABILITY = dplyr::case_when(
      ACTIONABILITY_TIER == 1 ~ "Strong significance",
      ACTIONABILITY_TIER == 2 ~ "Potential significance",
      ACTIONABILITY_TIER == 3 ~ "Uncertain significance",
      TRUE ~ as.character(NA)
    ))

  if(NROW(biomarker_items) > 0){
    results_amp_asco[['biomarker_evidence']][['items_all']] <-
      biomarker_items
  }

  if(NROW(biomarker_items_hires) > 0){
    results_amp_asco[['biomarker_evidence']][['items']] <-
      biomarker_items_hires |>
      dplyr::rename(ACTIONABILITY_TIER = .data$AMP_ASCO_TIER) |>
      dplyr::mutate(ACTIONABILITY_FRAMEWORK = "AMP_ASCO_CAP") |>
      dplyr::mutate(ACTIONABILITY = dplyr::case_when(
        ACTIONABILITY_TIER == 1 ~ "Strong significance",
        ACTIONABILITY_TIER == 2 ~ "Potential significance",
        ACTIONABILITY_TIER == 3 ~ "Uncertain significance",
        TRUE ~ as.character(NA)
      ))
  }

  if(NROW(tier_classification) > 0){
    results_amp_asco[['biomarker_evidence']][['tier_classification']] <-
      tier_classification |>
      dplyr::rename(ACTIONABILITY_TIER = .data$AMP_ASCO_TIER) |>
      dplyr::mutate(ACTIONABILITY_FRAMEWORK = "AMP_ASCO_CAP") |>
      dplyr::mutate(ACTIONABILITY = dplyr::case_when(
        ACTIONABILITY_TIER == 1 ~ "Strong significance",
        ACTIONABILITY_TIER == 2 ~ "Potential significance",
        ACTIONABILITY_TIER == 3 ~ "Uncertain significance",
        TRUE ~ as.character(NA)
      ))
  }

  return(results_amp_asco)

}
