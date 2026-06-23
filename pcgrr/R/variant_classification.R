#' AMP/ASCO/CAP tier classification for somatic variants in cancer
#'
#' Function that assigns tier classifications to genes subject to somatic CNA
#' (amplifications and deletions), somatic SNVs/InDels, or RNA fusions, based on
#' the presence and strength of biomarker evidence associated with records in
#' the variant set. Tier classification is based on the AMP/ASCO/CAP guidelines
#' for somatic variant interpretation in cancer (Li et al., J Mol Diagn. 2017).
#' The function also considers variant properties associated with oncogenes and
#' tumor suppressor genes (e.g. low MAF, coding status), which are used to
#' assign tier 3 status to variants with uncertain clinical significance.
#'
#' @param vartype variant type ('snv_indel', 'cna', 'fusion')
#' @param clinical_significance character indicate the clinical significance
#' types of evidence used for tiering of predictive evidence. Possible values:
#' "therapeutic_sensitivity", "therapeutic_resistance", "prognostic_poor",
#' "prognostic_better", "diagnostic_positive", "diagnostic_negative"
#' @param primary_site primary tumor site
#' @param biomarker_mapping_confidence confidence level of variant-biomarker
#' mapping resolution  (e.g. 'high' or 'medium'). 'High' indicates matches
#' at genomic/hgvsp/hgvsc level, 'medium' indicates matches at gene/exon level
#' Note that that matches at gene/exon level also consider variant properties
#' (loss-of-function, hotspot overlap etc) to avoid noise from presumably
#' less impactful variants that are likely not relevant for the biomarker
#' relationship
#' @param var_df data frame with variants (SNVs/InDels, CNAs, fusions)
#' @param biomarker_items data frame with biomarker evidence items
#'
#' @export
assign_amp_asco_cap_tiers <- function(
    vartype = "snv_indel",
    clinical_significance = "therapeutic_sensitivity",
    primary_site = "Any",
    biomarker_mapping_confidence = "medium",
    var_df = NULL,
    biomarker_items = NULL) {

  invisible(assertthat::assert_that(
    !is.null(var_df) & is.data.frame(var_df),
    msg = paste0("Argument var_df needs be of type data.frame")))
  invisible(
    assertthat::assert_that(
      vartype %in% c("snv_indel", "cna", "fusion"),
      msg = paste0("Argument 'vartype' needs to be one of
                   'snv_indel', 'cna' or 'fusion'"))
  )

  ## check that clinical_significance is a character of length 1
  invisible(assertthat::assert_that(
    is.character(clinical_significance) &
      length(clinical_significance) == 1,
    msg = "Argument 'clinical_significance' needs to be a character of length 1")
  )

  invisible(assertthat::assert_that(
    clinical_significance %in% names(bm_categories),
    msg = paste0("Argument 'clinical_significance' needs to be one of: ",
                 paste0(names(bm_categories), collapse = ", '")))
  )

  etype_for_tiering <-
    bm_categories[[clinical_significance]]$etype

  invisible(assertthat::assert_that(
    primary_site %in% tumor_sites,
    msg = paste0("Argument 'primary_site' needs to be one of: ",
                 paste0(tumor_sites, collapse = ", ")))
  )

  clnsig_categories_for_tiering <-
    bm_categories[[clinical_significance]]$clnsig

  if (!is.null(biomarker_items)) {

    invisible(assertthat::assert_that(
      is.data.frame(biomarker_items),
      msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

    if (NROW(biomarker_items) > 0) {

      invisible(assertable::assert_colnames(
        biomarker_items,
        c("BM_CLINICAL_SIGNIFICANCE",
          "BM_EVIDENCE_DIRECTION",
          "BM_EVIDENCE_TYPE",
          "BM_MAPPING_CONFIDENCE"),
        only_colnames = FALSE, quiet = T
      ))

      ## only consider biomarker evidence items with specified
      ## clinical significance
      biomarker_items <- biomarker_items |>
        dplyr::filter(
          .data$BM_CLINICAL_SIGNIFICANCE %in%
            clnsig_categories_for_tiering &

            ## only consider positive evidence items for
            ## tier classification (i.e. do not consider
            ## evidence items with 'Does Not Support' evidence
            ## direction)
            (.data$BM_EVIDENCE_DIRECTION == "Supports" |
            is.na(.data$BM_EVIDENCE_DIRECTION)))

      if (NROW(biomarker_items) > 0) {
        ## do not consider biomarkers with confidence below
        ## specified confidence level for tier classification
        if (biomarker_mapping_confidence == "medium") {
          biomarker_items <- biomarker_items |>
            dplyr::filter(
              .data$BM_MAPPING_CONFIDENCE %in% c("high", "medium")
            )
        }else{
          if (biomarker_mapping_confidence == "high") {
            biomarker_items <- biomarker_items |>
              dplyr::filter(
                .data$BM_MAPPING_CONFIDENCE == "high"
              )
          }
        }

        ## for fusion variants, do not consider biomarker evidence items with
        ## 'cgi' as source database, as these are poorly mapped, only
        ## at the level of individual genes, no 5' and 3' gene partner information'
        if(vartype == "fusion"){
          biomarker_items <- biomarker_items |>
            dplyr::filter(
              .data$BM_SOURCE_DB != "cgi"
            )
        }

      }
    }
  }

  if (NROW(biomarker_items) == 0) {
    log4r_debug(
      paste0("No biomarker evidence items with specified confidence ",
      "level found for tier classification"))
  }


  if (vartype == "snv_indel") {

    variants_tier_classified <- assign_variant_tiers_snv_indel(
      biomarker_items = biomarker_items,
      biomarker_mapping_confidence = biomarker_mapping_confidence,
      var_df = var_df,
      etype_for_tiering = etype_for_tiering,
      primary_site = primary_site
    )
  }else{
    if (vartype == "cna") {
      variants_tier_classified <- assign_variant_tiers_cna(
        biomarker_items = biomarker_items,
        biomarker_mapping_confidence = biomarker_mapping_confidence,
        var_df = var_df,
        etype_for_tiering = etype_for_tiering,
        primary_site = primary_site
      )
    }else{
      if (vartype == "fusion") {
        variants_tier_classified <- assign_variant_tiers_fusion(
          biomarker_items = biomarker_items,
          biomarker_mapping_confidence = biomarker_mapping_confidence,
          var_df = var_df,
          etype_for_tiering = etype_for_tiering,
          primary_site = primary_site
        )
      }
    }
  }

  ## debug the number of variants in 'variants_tier_classified'
  ## and the number of variants per unique tier in 'ACTIONABILITY_TIER'

  ## first log clinical_significance and variant type for context
  log4r_debug(paste0(
    "assign_amp_asco_cap_tiers - clinical_significance: ",
    clinical_significance,
    "; vartype: ", vartype
  ))
  log4r_debug(paste0(
    "assign_amp_asco_cap_tiers - number of variants in variants_tier_classified: ",
    NROW(variants_tier_classified)))
   log4r_debug(paste0(
     "assign_amp_asco_cap_tiers - number of variants per ACTIONABILITY_TIER: ",
     paste0(
       variants_tier_classified |>
         dplyr::group_by(.data$ACTIONABILITY_TIER) |>
         dplyr::summarise(n = dplyr::n()) |>
         dplyr::arrange(.data$ACTIONABILITY_TIER) |>
         dplyr::mutate(ACTIONABILITY_TIER = as.character(.data$ACTIONABILITY_TIER)),
       collapse = ", "
     )))


  if (primary_site != "Any") {
    biomarker_items <- assign_bm_tier_support_ttspecific(
      variants_tier_classified = variants_tier_classified,
      biomarker_items = biomarker_items,
      vartype = vartype,
      primary_site = primary_site,
      etype_for_tiering = etype_for_tiering
    )
  }else{
    biomarker_items <- assign_bm_tier_support_ttagnostic(
      variants_tier_classified = variants_tier_classified,
      biomarker_items = biomarker_items,
      vartype = vartype,
      etype_for_tiering = etype_for_tiering
    )
  }

  results_amp_asco <- list()
  results_amp_asco[['bm_evidence']] <- list()
  #results_amp_asco[[]]
  results_amp_asco[['variant']] <- var_df |>
    dplyr::left_join(
      variants_tier_classified,
      by = c("VAR_ID",
             "VARIANT_CLASS",
             "ENTREZGENE")
    )

  variants_tier_classified$ACTIONABILITY <- NULL

  ## if primary_site != "Any" and clinical significance not related to
  ## drug sensitivity or resistance, then limit biomarker evidence items to
  ## those with matching tumor-type hits
  if((clinical_significance == "prognostic_better" |
      clinical_significance == "prognostic_poor" |
      clinical_significance == "therapeutic_resistance" |
      clinical_significance == "diagnostic_positive" |
      clinical_significance == "diagnostic_negative") &
     primary_site != "Any" &
     NROW(biomarker_items) > 0) {
    biomarker_items <- biomarker_items |>
      dplyr::filter(.data$BM_PRIMARY_SITE == primary_site)

    ## log how many biomarker items / clinical significance items for the particular site
    # log4r_info(paste0(
    #   "assign_amp_asco_cap_tiers - ", clinical_significance,
    #   " - number of biomarker items with primary site match: ",
    #   NROW(biomarker_items)))
    #
    # ## log how many biomarker items / clinical significance items for the particular site
    # log4r_info(paste0(
    #   "assign_amp_asco_cap_tiers - ", clinical_significance,
    #   " - number of variants tier classified: ",
    #   NROW(variants_tier_classified)))

    if(NROW(biomarker_items) > 0 &
       NROW(variants_tier_classified) > 0){
      variants_tier_classified <- variants_tier_classified |>
        dplyr::filter(
          .data$VAR_ID %in% biomarker_items$VAR_ID &
            .data$ENTREZGENE %in% biomarker_items$ENTREZGENE &
            .data$VARIANT_CLASS %in% biomarker_items$VARIANT_CLASS
        )

      ## log remaining rows in variants_tier_classified after filtering for primary site match
      # log4r_info(paste0(
      #   "assign_amp_asco_cap_tiers - ", clinical_significance,
      #   " - number of variants in variants_tier_classified after filtering for primary site match: ",
      #   NROW(variants_tier_classified)))
    }

  }

  if(NROW(biomarker_items) == 0){
    results_amp_asco[['bm_evidence']][['classification']] <- data.frame()
  }else{
    results_amp_asco[['bm_evidence']][['classification']] <-
      variants_tier_classified |>
      dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
      dplyr::distinct()
  }
  results_amp_asco[['bm_evidence']][['eitems']] <-
    biomarker_items

  return(results_amp_asco)

}

#' Assign tiers of clinical significance (AMP/ASCO/CAP framework) to
#' somatic CNAs
#'
#' Function that assigns tiers of clinical significance (AMP/ASCO/CAP
#' framework) to somatic CNAs based on biomarker evidence items.
#' The function considers the strength of evidence (evidence levels)
#' and the match between biomarker site and primary site of query tumor.
#' The function also considers gene properties - oncogenes
#' and tumor suppressor genes - which are used to assign tier 3 status to
#' variants with uncertain clinical significance.
#'
#' @param primary_site primary tumor site of query tumor
#' (e.g. 'Lung', 'Breast', 'Any')
#' @param biomarker_mapping_confidence confidence level of variant-biomarker mapping
#' resolution  (e.g. 'high' or 'medium')
#' @param var_df data frame with somatic CNAs
#' @param etype_for_tiering evidence type(s) to consider for tier classification
#' (e.g. 'predictive', 'prognostic', 'diagnostic')
#' @param biomarker_items data frame with biomarker evidence items
#' @return data frame with tier classifications for somatic CNAs
#' based on biomarker evidence items and variant properties
#'
#' @export
#'
assign_variant_tiers_cna <- function(
    primary_site = "Any",
    biomarker_mapping_confidence = "medium",
    var_df = NULL,
    etype_for_tiering = c("predictive"),
    biomarker_items = NULL) {

  invisible(assertthat::assert_that(
    !is.null(var_df) &
      is.data.frame(var_df),
    msg = paste0(
      "Argument var_df needs be of type data.frame")))

  invisible(assertable::assert_colnames(
    var_df, c("VAR_ID",
              "VARIANT_CLASS",
              "ENTREZGENE",
              "TUMOR_SUPPRESSOR",
              "ONCOGENE"),
    only_colnames = FALSE, quiet = TRUE))

  invisible(assertthat::assert_that(
    primary_site %in% tumor_sites,
    msg = paste0(
      "Argument 'primary_site' needs to be one of: ",
      paste0(tumor_sites, collapse = ", ")))
  )

  invisible(assertthat::assert_that(
    !is.null(biomarker_items) & is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

  log4r_debug(paste0(
    "assign_variant_tiers_cna - etype_for_tiering: ",
    paste(etype_for_tiering, collapse = ", ")))

  variants_tier_classified <- data.frame()

  if (primary_site != "Any") {
    variants_tier_classified <- assign_variant_top_tiers_ttspecific(
      biomarker_items = biomarker_items,
      biomarker_mapping_confidence = biomarker_mapping_confidence,
      var_df = var_df,
      etype_for_tiering = etype_for_tiering,
      primary_site = primary_site
    )
  }else{
    variants_tier_classified <- assign_variant_top_tiers_ttagnostic(
      biomarker_items = biomarker_items,
      biomarker_mapping_confidence = biomarker_mapping_confidence,
      var_df = var_df,
      etype_for_tiering = etype_for_tiering
    )
  }

  ## log length gene variant records (VAR_ID, ENTREZGENE, VARIANT_CLASS)
  log4r_debug(paste0(
    "assign_variant_tiers_cna - number of variants for tier classification: ",
    NROW(variants_tier_classified)))


  other_variant_properties <- var_df |>
    dplyr::select(
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "TUMOR_SUPPRESSOR",
        "ONCOGENE")
    )

  invisible(assertable::assert_colnames(
    variants_tier_classified,
    c("VAR_ID",
      "VARIANT_CLASS",
      "ENTREZGENE",
      "ACTIONABILITY_TIER"),
    only_colnames = TRUE, quiet = TRUE))

  ## Assign low tier status to variants with uncertain clinical significance
  ## based on gene properties
  ## 1) amplified oncogenes, OR
  ## 2) lost tumor suppressor genes (het or homozygous deletions)

  variants_tier_classified <- variants_tier_classified |>
    dplyr::left_join(
      other_variant_properties,
      by = c("VAR_ID",
             "VARIANT_CLASS",
             "ENTREZGENE")) |>

    dplyr::mutate(
      ACTIONABILITY_TIER = dplyr::if_else(
        (!is.na(.data$TUMOR_SUPPRESSOR) &
           .data$TUMOR_SUPPRESSOR == TRUE &
           (.data$VARIANT_CLASS == "homdel" |
              .data$VARIANT_CLASS == "hemdel" |
              .data$VARIANT_CLASS == "hetdel")) |
          (!is.na(.data$ONCOGENE) &
             .data$ONCOGENE == TRUE &
             .data$VARIANT_CLASS == "amplification") &
          is.na(.data$ACTIONABILITY_TIER),
        as.integer(3),
        as.integer(.data$ACTIONABILITY_TIER))) |>
    dplyr::arrange(.data$ACTIONABILITY_TIER) |>
    dplyr::distinct() |>
    dplyr::mutate(
      ACTIONABILITY = dplyr::case_when(
        ACTIONABILITY_TIER == 1 ~ "Strong significance",
        ACTIONABILITY_TIER == 2 ~ "Potential significance",
        ACTIONABILITY_TIER == 3 ~ "Uncertain significance",
        TRUE ~ as.character(NA)
      )
    ) |>
    dplyr::select(
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "ACTIONABILITY_TIER",
        "ACTIONABILITY")
    )

  return(variants_tier_classified)

}

#' Assign tiers of clinical significance (AMP/ASCO/CAP framework) to
#' RNA fusions
#'
#' Function that assigns tiers of clinical significance (AMP/ASCO/CAP
#' framework) to RNA fusions based on biomarker evidence items.
#' The function considers the strength of evidence (evidence levels)
#' and the match between biomarker site and primary site of query tumor.
#' The function also considers oncogene properties of fusion partners to
#' assign tier 3 status to variants with uncertain clinical significance.
#'
#' @param primary_site primary tumor site of query tumor
#' (e.g. 'Lung', 'Breast', 'Any')
#' @param biomarker_mapping_confidence confidence level of variant-biomarker mapping
#' resolution  (e.g. 'high' or 'medium')
#' @param var_df data frame with somatic RNA fusions
#' @param etype_for_tiering evidence type(s) to consider for tier classification
#' (e.g. 'predictive', 'prognostic', 'diagnostic')
#' @param biomarker_items data frame with biomarker evidence items
#' @return data frame with tier classifications for somatic RNA fusions
#' based on biomarker evidence items and variant properties
#'
#' @export
#'
assign_variant_tiers_fusion <- function(
    primary_site = "Any",
    biomarker_mapping_confidence = "medium",
    var_df = NULL,
    etype_for_tiering = c("predictive"),
    biomarker_items = NULL) {

  invisible(assertthat::assert_that(
    !is.null(var_df) & is.data.frame(var_df),
    msg = paste0("Argument var_df needs be of type data.frame")))

  invisible(assertable::assert_colnames(
    var_df, c("VAR_ID",
              "VARIANT_CLASS",
              "ENTREZGENE",
              "ONCOGENE_3P",
              "ONCOGENE_5P"),
    only_colnames = FALSE, quiet = TRUE))

  invisible(assertthat::assert_that(
    primary_site %in% tumor_sites,
    msg = paste0("Argument 'primary_site' needs to be one of: ",
                 paste0(tumor_sites, collapse = ", ")))
  )

  invisible(assertthat::assert_that(
    !is.null(biomarker_items) & is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

  log4r_debug(paste0(
    "assign_variant_tiers_fusion - etype_for_tiering: ",
    paste(etype_for_tiering, collapse = ", ")))

  variants_tier_classified <- data.frame()

  if (primary_site != "Any") {
    variants_tier_classified <- assign_variant_top_tiers_ttspecific(
      biomarker_items = biomarker_items,
      biomarker_mapping_confidence = biomarker_mapping_confidence,
      var_df = var_df,
      etype_for_tiering = etype_for_tiering,
      primary_site = primary_site
    )
  }else{
    variants_tier_classified <- assign_variant_top_tiers_ttagnostic(
      biomarker_items = biomarker_items,
      biomarker_mapping_confidence = biomarker_mapping_confidence,
      var_df = var_df,
      etype_for_tiering = etype_for_tiering
    )
  }

  ## log length gene variant records (VAR_ID, ENTREZGENE, VARIANT_CLASS)
  log4r_debug(paste0(
    "assign_variant_tiers_fusion - number of variants for tier classification: ",
    NROW(variants_tier_classified)))

  other_variant_properties <- var_df |>
    dplyr::select(
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "MITDB_NUM_EVIDENCE",
        "ONCOGENE_3P",
        "ONCOGENE_5P")
        )

  invisible(assertable::assert_colnames(
    variants_tier_classified,
    c("VAR_ID",
      "VARIANT_CLASS",
      "ENTREZGENE",
      "ACTIONABILITY_TIER"),
    only_colnames = TRUE, quiet = TRUE))

  ## Assign low tier status to variants with uncertain clinical significance
  ## based oncogene status of fusion partners

  variants_tier_classified <- variants_tier_classified |>
    dplyr::left_join(
      other_variant_properties,
      by = c("VAR_ID",
             "VARIANT_CLASS",
             "ENTREZGENE")) |>
    dplyr::mutate(ACTIONABILITY_TIER = dplyr::if_else(
      is.na(.data$ACTIONABILITY_TIER) &
        ((!is.na(.data$ONCOGENE_3P) &
            .data$ONCOGENE_3P == TRUE) |
           (!is.na(.data$ONCOGENE_5P) &
              .data$ONCOGENE_5P == TRUE) &
           (!is.na(MITDB_NUM_EVIDENCE) &
              .data$MITDB_NUM_EVIDENCE > 0)
        ),
      as.integer(3),
      as.integer(.data$ACTIONABILITY_TIER)
    )) |>
    dplyr::mutate(ACTIONABILITY_TIER = dplyr::if_else(
      is.na(.data$ACTIONABILITY_TIER),
      as.integer(5),
      as.integer(.data$ACTIONABILITY_TIER))) |>
    dplyr::arrange(.data$ACTIONABILITY_TIER) |>
    dplyr::distinct() |>
    dplyr::mutate(
      ACTIONABILITY = dplyr::case_when(
        ACTIONABILITY_TIER == 1 ~ "Strong significance",
        ACTIONABILITY_TIER == 2 ~ "Potential significance",
        ACTIONABILITY_TIER == 3 ~ "Uncertain significance",
        TRUE ~ as.character(NA)
      )
    ) |>
    dplyr::select(
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "ACTIONABILITY_TIER",
        "ACTIONABILITY")
    )

  return(variants_tier_classified)

}

#' Assign tiers of clinical significance (AMP/ASCO/CAP framework) to
#' somatic SNVs/InDels
#'
#' Function that assigns tiers of clinical significance (AMP/ASCO/CAP
#' framework) to somatic SNVs/InDels based on biomarker evidence items.
#' The function considers the strength of evidence (evidence levels)
#' and the match between biomarker site and primary site of query tumor.
#' The function also considers variant properties associated with oncogenes
#' and tumor suppressor genes (e.g. low MAF, coding status), which are used
#' to assign tier 3 status to variants with uncertain clinical significance.
#'
#' @param primary_site primary tumor site of query tumor
#' (e.g. 'Lung', 'Breast', 'Any')
#' @param biomarker_mapping_confidence confidence level of variant-biomarker mapping
#' resolution  (e.g. 'high' or 'medium')
#' @param var_df data frame with somatic SNV/InDel variants
#' @param etype_for_tiering evidence type(s) to consider for tier classification
#' (e.g. 'predictive', 'prognostic', 'diagnostic')
#' @param biomarker_items data frame with biomarker evidence items
#' @return data frame with tier classifications for somatic SNVs/InDels
#' based on biomarker evidence items and variant properties
#'
#' @export
#'
assign_variant_tiers_snv_indel <- function(
    primary_site = "Any",
    biomarker_mapping_confidence = "medium",
    var_df = NULL,
    etype_for_tiering = c("predictive"),
    biomarker_items = NULL) {

  invisible(assertthat::assert_that(
    !is.null(var_df) & is.data.frame(var_df),
    msg = paste0("Argument var_df needs be of type data.frame")))

  invisible(assertable::assert_colnames(
    var_df, c("TUMOR_SUPPRESSOR",
              "VAR_ID",
              "VARIANT_CLASS",
              "ENTREZGENE",
              "ONCOGENE",
              "ONCOGENICITY",
              "gnomADe_AF",
              "CODING_STATUS"),
    only_colnames = FALSE, quiet = TRUE))

  invisible(assertthat::assert_that(
    primary_site %in% tumor_sites,
    msg = paste0("Argument 'primary_site' needs to be one of: ",
                 paste0(tumor_sites, collapse = ", ")))
  )

  invisible(assertthat::assert_that(
    is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

  log4r_debug(paste0(
    "assign_variant_tiers_snv_indel - etype_for_tiering: ",
    paste(etype_for_tiering, collapse = ", ")))

  ## log length gene variant records (VAR_ID, ENTREZGENE, VARIANT_CLASS)
  log4r_debug(paste0(
    "assign_variant_tiers_snv_indel - number of variants for tier classification: ",
    NROW(var_df)))

  variants_tier_classified <- data.frame()

  if (primary_site != "Any") {
    variants_tier_classified <- assign_variant_top_tiers_ttspecific(
      biomarker_items = biomarker_items,
      biomarker_mapping_confidence = biomarker_mapping_confidence,
      var_df = var_df,
      etype_for_tiering = etype_for_tiering,
      primary_site = primary_site
    )
  }else{
    variants_tier_classified <- assign_variant_top_tiers_ttagnostic(
      biomarker_items = biomarker_items,
      biomarker_mapping_confidence = biomarker_mapping_confidence,
      var_df = var_df,
      etype_for_tiering = etype_for_tiering
    )
  }

  ## log length gene variant records (VAR_ID, ENTREZGENE, VARIANT_CLASS)
  log4r_debug(paste0(
    "assign_variant_tiers_snv_indel - number of variants for tier classification: ",
    NROW(variants_tier_classified)))

  other_variant_properties <- var_df |>
    dplyr::select(
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "TUMOR_SUPPRESSOR",
        "ONCOGENE",
        "ONCOGENICITY",
        "gnomADe_AF",
        "CODING_STATUS")
    )

  invisible(assertable::assert_colnames(
    variants_tier_classified,
    c("VAR_ID",
      "VARIANT_CLASS",
      "ENTREZGENE",
      "ACTIONABILITY_TIER"),
    only_colnames = TRUE, quiet = TRUE))

  ## Assign low tier status to variants with uncertain clinical significance
  ## based on variant properties associated with oncogenes and tumor
  ## suppressor genes (e.g. low MAF, coding status)

  variants_tier_classified <- variants_tier_classified |>
    dplyr::left_join(
      other_variant_properties,
      by = c("VAR_ID",
             "VARIANT_CLASS",
             "ENTREZGENE")) |>
    dplyr::distinct() |>

    ## Tumor suppressor/oncogene mutations (tier 3)
    ## - low MAF (< 0.001)
    dplyr::mutate(ACTIONABILITY_TIER = dplyr::if_else(
      ((!is.na(.data$TUMOR_SUPPRESSOR) &
          .data$TUMOR_SUPPRESSOR == TRUE) |
         (!is.na(.data$ONCOGENE) &
            .data$ONCOGENE == TRUE)) &
        (is.na(.data$gnomADe_AF) |
           .data$gnomADe_AF < 0.001) &
        .data$CODING_STATUS == "coding" &
        (.data$ONCOGENICITY != "Benign" |
           .data$ONCOGENICITY != "Likely Benign") &
        is.na(.data$ACTIONABILITY_TIER),
      as.integer(3),
      as.integer(.data$ACTIONABILITY_TIER)
    )) |>

    ## Tier 4 - as specified in AMP/ASCO guidelines are
    ## not considered here

    ## Tier 5 - other coding mutations
    dplyr::mutate(ACTIONABILITY_TIER = dplyr::if_else(
      is.na(.data$ACTIONABILITY_TIER) &
        .data$CODING_STATUS == "coding",
      as.integer(5),
      as.integer(.data$ACTIONABILITY_TIER)
    )) |>

    ## Tier 6 - other non-coding mutations
    dplyr::mutate(ACTIONABILITY_TIER = dplyr::if_else(
      is.na(.data$ACTIONABILITY_TIER) &
        .data$CODING_STATUS == "noncoding",
      as.integer(6),
      as.integer(.data$ACTIONABILITY_TIER)
    )) |>
    dplyr::arrange(.data$ACTIONABILITY_TIER) |>
    dplyr::distinct() |>
    dplyr::arrange(
      .data$ACTIONABILITY_TIER
    ) |>
    dplyr::mutate(
      ACTIONABILITY = dplyr::case_when(
        ACTIONABILITY_TIER == 1 ~ "Strong significance",
        ACTIONABILITY_TIER == 2 ~ "Potential significance",
        ACTIONABILITY_TIER == 3 ~ "Uncertain significance",
        TRUE ~ as.character(NA)
      )
    ) |>
    dplyr::select(
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "ACTIONABILITY_TIER",
        "ACTIONABILITY")
    )

  return(variants_tier_classified)

}


#' Assign top tiers of clinical significance (AMP/ASCO/CAP framework) to
#' variants based on data from biomarker evidence items (tumor-type agnostic query)
#'
#' The function considers the strength of evidence (evidence levels) and
#' the match between biomarker site and primary site of query tumor.
#'
#' @param biomarker_items data frame with biomarker evidence items
#' @param biomarker_mapping_confidence confidence level of variant-biomarker mapping
#' resolution  (e.g. 'high' or 'medium')
#' @param var_df data frame with variants (SNVs/InDels, CNAs, fusions)
#' @param etype_for_tiering evidence type(s) to consider for tier classification
#' (e.g. 'predictive', 'prognostic', 'diagnostic')
#' @return data frame with top tier classifications for variants based on biomarker evidence items
#'
#' @export
#'
assign_variant_top_tiers_ttagnostic <- function(
    biomarker_items = NULL,
    biomarker_mapping_confidence = "high",
    var_df = NULL,
    etype_for_tiering = c("predictive")) {

  valid_etypes <-
    bm_evidence$types[
      !bm_evidence$types %in% c("oncogenic","functional")]

  invisible(assertthat::assert_that(
    is.character(etype_for_tiering),
    msg = "etype_for_tiering must be a character vector"))

  invisible(assertthat::assert_that(
    length(etype_for_tiering) > 0,
    msg = "etype_for_tiering must contain at least one value"))

  invisible(assertthat::assert_that(
    all(etype_for_tiering %in% valid_etypes),
    msg = sprintf("etype_for_tiering must be one or more of: %s",
                  paste(valid_etypes, collapse = ", "))))

  assertthat::assert_that(!any(duplicated(etype_for_tiering)),
                          msg = "etype_for_tiering contains duplicate values")

  invisible(assertthat::assert_that(
    biomarker_mapping_confidence %in% c("high", "medium"),
    msg = paste0("Argument 'biomarker_mapping_confidence' needs ",
                 "to be one of 'high' or 'medium'"))
  )

  invisible(assertthat::assert_that(
    is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

  invisible(assertthat::assert_that(
    is.data.frame(var_df) & !is.null(var_df),
    msg = paste0("Argument 'var_df' needs be of type data.frame")))

  invisible(assertable::assert_colnames(
    var_df, c("VAR_ID",
              "ENTREZGENE",
              "VARIANT_CLASS"),
    only_colnames = FALSE, quiet = TRUE))

  ## log length gene variant records (VAR_ID, ENTREZGENE, VARIANT_CLASS)
  log4r_debug(paste0(
    "assign_variant_top_tiers_ttagnostic - number of variants for tier classification: ",
    NROW(var_df)))

  log4r_debug(paste0(
    "assign_variant_top_tiers_ttagnostic - etype_for_tiering: ",
    paste(etype_for_tiering, collapse = ", ")))

  variants_tier_classified <- data.frame()

  ## Assign tiers based on biomarker evidence items
  ## (if any biomarker items with specified confidence level are available)
  if (!is.null(biomarker_items) & NROW(biomarker_items) > 0) {

    invisible(assertable::assert_colnames(
      biomarker_items,
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "BM_EVIDENCE_LEVEL",
        "BM_EVIDENCE_TYPE",
        "BM_PRIMARY_SITE"),
      only_colnames = FALSE, quiet = TRUE))

    ## For fusions, the ENTREZGENE in biomarker_items is the single-gene DB
    ## value (e.g. "238" for ALK fusions) while var_df carries the resolved
    ## pair (e.g. "238::238"). Grouping and joining on ENTREZGENE silently
    ## fails for single-gene fusion biomarkers, leaving every fusion with NA
    ## tier. VAR_ID uniquely identifies each fusion event in the sample, so
    ## use (VAR_ID, VARIANT_CLASS) only as the grouping/join key for fusions.
    is_fusion <- isTRUE("fusion" %in% unique(var_df$VARIANT_CLASS))
    tier_group_cols <- if (is_fusion) {
      c("VAR_ID", "VARIANT_CLASS")
    } else {
      c("VAR_ID", "VARIANT_CLASS", "ENTREZGENE")
    }

    biomarkers_tier_classified <-
      biomarker_items |>
      dplyr::select(
        c("VAR_ID",
          "VARIANT_CLASS",
          "ENTREZGENE",
          "BM_EVIDENCE_LEVEL",
          "BM_EVIDENCE_TYPE",
          "BM_PRIMARY_SITE")) |>
      dplyr::distinct() |>
      dplyr::mutate(ACTIONABILITY_TIER = dplyr::case_when(

        ## A) Biomarker site is pan-cancer ('Any')
        ## B) strong evidence - evidence levels A/B
        ## --> TIER 1
        .data$BM_PRIMARY_SITE == 'Any' &
          tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL,
            bm_evidence$strong_regex
          ) ~ as.integer(1),

        ## A) Biomarker site is _NOT_ pan-cancer ('Any')
        ## B) strong evidence - evidence levels A/B
        ## --> TIER 2

        .data$BM_PRIMARY_SITE != "Any" &
          tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL,
            bm_evidence$strong_regex
          ) ~ as.integer(2),

        ## A) weak evidence - evidence levels (C/D/E)
        ## B) Pan-cancer biomarker or tumor-type specific biomarker
        ## --> TIER 3 (any weak evidence is considered tier 3,
        ##     regardless of tumor-type specificity)
        tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL,
            bm_evidence$weak_regex
          ) ~ as.integer(3),
        TRUE ~ as.integer(100)

      )) |>

      ## Get top tier level for any given variant.
      ## For fusions group on (VAR_ID, VARIANT_CLASS) only — see note above.
      dplyr::group_by(dplyr::across(dplyr::all_of(tier_group_cols))) |>
      dplyr::summarise(
        ACTIONABILITY_TIER = min(
          .data$ACTIONABILITY_TIER, na.rm = TRUE),
        .groups = "drop") |>
      dplyr::mutate(ACTIONABILITY_TIER = dplyr::if_else(
        .data$ACTIONABILITY_TIER == 100,
        as.integer(NA),
        as.integer(.data$ACTIONABILITY_TIER)
      ))

    if (NROW(biomarkers_tier_classified) > 0 &
       NROW(var_df) > 0) {

      variants_tier_classified <- var_df |>
        dplyr::select(
          c("VAR_ID",
            "VARIANT_CLASS",
            "ENTREZGENE")
        ) |>
        dplyr::left_join(
          dplyr::select(
            biomarkers_tier_classified,
            dplyr::all_of(c(tier_group_cols, "ACTIONABILITY_TIER"))
          ),
          by = tier_group_cols) |>
        dplyr::distinct()|>
        dplyr::arrange(
          .data$ACTIONABILITY_TIER
        )
    }

  }else{
    ## if no biomarker evidence items with specified confidence level are
    ## available, assign NA tier status to all variants
    log4r_debug(
      paste0("No biomarker evidence items with specified confidence ",
             "level found for tier classification"))

    variants_tier_classified <- var_df |>
      dplyr::select(
        c("VAR_ID",
          "VARIANT_CLASS",
          "ENTREZGENE")
      ) |>
      dplyr::mutate(
        ACTIONABILITY_TIER = as.integer(NA),
      ) |>
      dplyr::distinct()
  }

  return(variants_tier_classified)


}

#' Assign top tiers of clinical significance (AMP/ASCO/CAP framework) to
#' variants based on data from biomarker evidence items (tumor type-specific query)
#'
#' The function considers the strength of evidence (evidence levels) and
#' the match between biomarker site and primary site of query tumor.
#'
#' @param biomarker_items data frame with biomarker evidence items
#' @param biomarker_mapping_confidence confidence level of variant-biomarker mapping
#' resolution  (e.g. 'high' or 'medium')
#' @param var_df data frame with variants (SNVs/InDels, CNAs, fusions)
#' @param etype_for_tiering evidence type(s) to consider for tier classification
#' (e.g. 'predictive', 'prognostic', 'diagnostic')
#' @param primary_site primary tumor site of query tumor (e.g. 'Lung', 'Breast')
#' @return data frame with top tier classifications for variants based on biomarker evidence items
#'
#' @export
#'
assign_variant_top_tiers_ttspecific <- function(
    biomarker_items = NULL,
    biomarker_mapping_confidence = "high",
    var_df = NULL,
    etype_for_tiering = c("predictive"),
    primary_site = "Lung") {

  invisible(assertthat::assert_that(
    primary_site %in% tumor_sites,
    msg = paste0("Argument 'primary_site' needs to be one of: ",
                 paste0(tumor_sites, collapse = ", ")))
  )

  ## check that primary_site is not 'Any' for tumor-specific tiering function
  invisible(assertthat::assert_that(
    primary_site != "Any",
    msg = "Argument 'primary_site' cannot be 'Any' for tumor-specific tiering
    function 'assign_variant_top_tiers_ttspecific'"))


  valid_etypes <-
    bm_evidence$types[
      !bm_evidence$types %in% c("oncogenic","functional")]

  invisible(assertthat::assert_that(
    is.character(etype_for_tiering),
    msg = "etype_for_tiering must be a character vector"))

  invisible(assertthat::assert_that(
    length(etype_for_tiering) > 0,
    msg = "etype_for_tiering must contain at least one value"))

  invisible(assertthat::assert_that(
    all(etype_for_tiering %in% valid_etypes),
    msg = sprintf("etype_for_tiering must be one or more of: %s",
                  paste(valid_etypes, collapse = ", "))))

  assertthat::assert_that(!any(duplicated(etype_for_tiering)),
              msg = "etype_for_tiering contains duplicate values")

  invisible(assertthat::assert_that(
    biomarker_mapping_confidence %in% c("high", "medium"),
    msg = paste0("Argument 'biomarker_mapping_confidence' needs ",
                 "to be one of 'high' or 'medium'"))
  )
  invisible(assertthat::assert_that(
    is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

  invisible(assertthat::assert_that(
    is.data.frame(var_df) & !is.null(var_df),
    msg = paste0("Argument 'var_df' needs be of type data.frame")))

  invisible(assertable::assert_colnames(
    var_df, c("VAR_ID",
              "ENTREZGENE",
              "VARIANT_CLASS"),
    only_colnames = FALSE, quiet = TRUE))

  ## log length gene variant records (VAR_ID, ENTREZGENE, VARIANT_CLASS)
  # log4r_debug(paste0(
  #   "assign_variant_top_tiers_ttspecific - number of variants for tier classification: ",
  #   NROW(var_df)))
  #
  # log4r_debug(paste0(
  #   "assign_variant_top_tiers_ttspecific - etype_for_tiering: ",
  #   paste(etype_for_tiering, collapse = ", ")))

  variants_tier_classified <- data.frame()

  ## Assign tiers based on biomarker evidence items
  ## (if any biomarker items with specified confidence level are available)
  if (!is.null(biomarker_items) & NROW(biomarker_items) > 0) {

    invisible(assertable::assert_colnames(
      biomarker_items,
      c("VAR_ID",
        "VARIANT_CLASS",
        "ENTREZGENE",
        "BM_EVIDENCE_LEVEL",
        "BM_EVIDENCE_TYPE",
        "BM_PRIMARY_SITE"),
      only_colnames = FALSE, quiet = TRUE))

    ## For fusions, the ENTREZGENE in biomarker_items is the single-gene DB
    ## value (e.g. "4916" for "NTRK3 fusions") while var_df carries the
    ## resolved pair (e.g. "123624::4916"). Grouping and joining on ENTREZGENE
    ## silently fails for single-gene fusion biomarkers, leaving every fusion
    ## with NA tier and falling through to the oncogene-partner Tier 3 fallback.
    ## VAR_ID uniquely identifies each fusion event, so use (VAR_ID, VARIANT_CLASS)
    ## only as the grouping/join key for fusions.
    is_fusion <- isTRUE(
      "fusion" %in% unique(var_df$VARIANT_CLASS))
    tier_group_cols <- if (is_fusion) {
      c("VAR_ID", "VARIANT_CLASS")
    } else {
      c("VAR_ID", "VARIANT_CLASS", "ENTREZGENE")
    }

    biomarkers_tier_classified <-
      biomarker_items |>
      dplyr::select(
        c("VAR_ID",
          "VARIANT_CLASS",
          "ENTREZGENE",
          "BM_EVIDENCE_LEVEL",
          "BM_EVIDENCE_TYPE",
          "BM_PRIMARY_SITE")) |>
      dplyr::distinct() |>
      dplyr::mutate(ACTIONABILITY_TIER = dplyr::case_when(

        ## A) Biomarker site _matches_ primary site of query tumor, OR
        ##.   pan-cancer biomarker (BM_PRIMARY_SITE = 'Any')
        ## B) strong evidence - evidence levels A/B
        ## --> TIER 1
        (.data$BM_PRIMARY_SITE == primary_site |
           .data$BM_PRIMARY_SITE == 'Any') &
          tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, bm_evidence$strong_regex
          ) ~ as.integer(1),

        ## A) Biomarker site _does not_ match primary site of query tumor
        ## B) Biomarker site is not 'Any' (non pan-cancer biomarker)
        ## C) strong evidence - evidence levels A/B
        ## --> TIER 2
        (.data$BM_PRIMARY_SITE != primary_site &
           .data$BM_PRIMARY_SITE != "Any") &
           tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
           stringr::str_detect(
             .data$BM_EVIDENCE_LEVEL, bm_evidence$strong_regex
           ) ~ as.integer(2),

        ## A) Biomarker site _matches_ primary site of query tumor
        ## B) weak evidence - evidence levels (C/D/E)
        ## --> TIER 2
        .data$BM_PRIMARY_SITE == primary_site &
          tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex
          ) ~ as.integer(2),

        ## A) Biomarker site _does not_ match primary site of query tumor
        ## B) Biomarker site is not pan-cancer ('Any')
        ## C) weak evidence - evidence levels (C/D/E)
        ## --> TIER 3
        .data$BM_PRIMARY_SITE != primary_site &
           .data$BM_PRIMARY_SITE != "Any" &
          tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex
          ) ~ as.integer(3),

        ## A) Pan-cancer biomarker (BM_PRIMARY_SITE = 'Any')
        ## B) Weak evidence - evidence levels (C/D/E)
        ## --> TIER 3 (pan-cancer weak is relevant signal, but at the
        ##     lowest tier — consistent with non-matching weak evidence)
        .data$BM_PRIMARY_SITE == "Any" &
          tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
          stringr::str_detect(
            .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex
          ) ~ as.integer(3),

        TRUE ~ as.integer(100)

      )) |>

      ## Get top tier level for any given variant.
      ## For fusions group on (VAR_ID, VARIANT_CLASS) only — see note above.
      dplyr::group_by(dplyr::across(dplyr::all_of(tier_group_cols))) |>
      dplyr::summarise(
        ACTIONABILITY_TIER = min(
          .data$ACTIONABILITY_TIER, na.rm = TRUE),
        .groups = "drop") |>
      dplyr::mutate(ACTIONABILITY_TIER = dplyr::if_else(
        .data$ACTIONABILITY_TIER == 100,
        as.integer(NA),
        as.integer(.data$ACTIONABILITY_TIER)
      ))

    if (NROW(biomarkers_tier_classified) > 0 &
       NROW(var_df) > 0) {

      ## log Number of rows in var_df
      log4r_debug(paste0(
        "assign_variant_top_tiers_ttspecific - number of variants in var_df: ",
        NROW(var_df)))

      variants_tier_classified <- var_df |>
        dplyr::select(
          c("VAR_ID",
            "VARIANT_CLASS",
            "ENTREZGENE")
        ) |>
        dplyr::left_join(
          dplyr::select(
            biomarkers_tier_classified,
            dplyr::all_of(c(tier_group_cols, "ACTIONABILITY_TIER"))
          ),
          by = tier_group_cols) |>
        dplyr::distinct()|>
        dplyr::arrange(
          .data$ACTIONABILITY_TIER
        )
    }

  }else{
    ## if no biomarker evidence items with specified confidence level are
    ## available, assign NA tier status to all variants
    log4r_debug(
      paste0("No biomarker evidence items with specified confidence ",
             "level found for tier classification"))

    variants_tier_classified <- var_df |>
      dplyr::select(
        c("VAR_ID",
          "VARIANT_CLASS",
          "ENTREZGENE")
      ) |>
      dplyr::mutate(
        ACTIONABILITY_TIER = as.integer(NA)
      ) |>
      dplyr::distinct()
  }

  return(variants_tier_classified)

}

#' For all biomarker evidence items with specified confidence level,
#' assign whether these are tier-defining or providing additional
#' support - tumor-type specific query
#'
#' The function considers the strength of evidence (evidence levels) and
#' the match between biomarker site and primary site of query tumor.
#' Mark evidence considered tier-defining, while other biomarker evidence items
#' as additional support
#'
#' @param variants_tier_classified data frame with variant tier classifications
#' based on support from biomarker evidence item
#' @param vartype variant type (e.g. 'snv_indel', 'cna', 'fusion')
#' @param primary_site primary tumor site of query tumor (e.g. 'Lung',
#' 'Breast', 'Any')
#' @param etype_for_tiering evidence type(s) to consider for tier classification
#' @param biomarker_items data frame with biomarker evidence items
#' @return data frame with biomarker evidence items classified as tier-defining
#' or providing additional support
#'
#' @export
#'
assign_bm_tier_support_ttspecific <- function(
    variants_tier_classified = NULL,
    vartype = NULL,
    primary_site = "Lung",
    etype_for_tiering = c("predictive"),
    biomarker_items = NULL) {

  invisible(assertthat::assert_that(
    !is.null(variants_tier_classified) &
      is.data.frame(variants_tier_classified),
    msg = paste0("Argument var_df needs be of type data.frame")))
  invisible(
    assertthat::assert_that(
      !is.null(vartype) & is.character(vartype),
      vartype %in% c("snv_indel", "cna", "fusion"),
      msg = paste0("Argument 'vartype' needs to be one of
                   'snv_indel', 'cna' or 'fusion'"))
  )

  invisible(assertthat::assert_that(
    primary_site %in% tumor_sites,
    msg = paste0("Argument 'primary_site' needs to be one of: ",
                 paste0(tumor_sites, collapse = ", ")))
  )

  ## check that primary_site is not 'Any' for tumor-specific tiering function
  invisible(assertthat::assert_that(
    primary_site != "Any",
    msg = paste0(
      "Argument 'primary_site' cannot be 'Any' for tumor-specific ",
      "actionability support function 'assign_bm_tier_support_ttspecific'")))

  invisible(assertthat::assert_that(
    is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

  assertable::assert_colnames(
    variants_tier_classified,
    c("VAR_ID",
      "VARIANT_CLASS",
      "ACTIONABILITY",
      "ACTIONABILITY_TIER"),
    only_colnames = FALSE, quiet = TRUE)

  log4r_debug(paste0(
    "assign_bm_tier_support_ttspecific - etype_for_tiering: ",
    paste(etype_for_tiering, collapse = ", ")))

  if (NROW(biomarker_items) > 0) {

    assertable::assert_colnames(
      biomarker_items,
      c("VAR_ID",
        "ENTREZGENE",
        "BM_EVIDENCE_LEVEL",
        "BM_PRIMARY_SITE"),
      only_colnames = FALSE, quiet = TRUE)

    ## For fusions, the ENTREZGENE in biomarker_items reflects the biomarker
    ## DB entry (e.g. "4916" for "NTRK3 fusions") while variants_tier_classified
    ## carries the resolved gene pair (e.g. "123624::4916"). Joining on ENTREZGENE
    ## would always fail for single-gene fusion biomarkers, leaving ACTIONABILITY_TIER
    ## as NA and producing "TNA" in the TIER column. VAR_ID already uniquely
    ## identifies each fusion event in the sample, so join on (VAR_ID, VARIANT_CLASS)
    ## only for fusions.
    tier_join_cols <- if (vartype == "fusion") {
      c("VAR_ID", "VARIANT_CLASS")
    } else {
      c("VAR_ID", "VARIANT_CLASS", "ENTREZGENE")
    }
    tier_lookup_cols <- c(tier_join_cols, "ACTIONABILITY_TIER")

    biomarker_items <- as.data.frame(
      biomarker_items |>
        dplyr::left_join(
          dplyr::select(
            variants_tier_classified,
            dplyr::all_of(tier_lookup_cols)
          ),
          by = tier_join_cols) |>
        dplyr::distinct() |>
        dplyr::select(
          c("VAR_ID",
            "VARIANT_CLASS",
            "ENTREZGENE",
            "ACTIONABILITY_TIER"),
          dplyr::everything()
        ) |>

        ## Now go through the actionability tier for each evidence item
        ## and set them as tier-defining or supporting based on the strength
        ## of evidence and tumor type match
        dplyr::mutate(
          BM_ACTIONABILITY_SUPPORT = dplyr::case_when(

            ## Pan-cancer weak evidence on tier-1/2 variants: the tier was
            ## set by stronger evidence — pan-cancer weak is additional context
            ## only. On tier-3 variants it falls through to the tier-3 rule
            ## below where it is correctly labelled tier-defining.
            .data$BM_PRIMARY_SITE == "Any" &
              .data$ACTIONABILITY_TIER %in% c(1L, 2L) &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex) ~ "additional",

            ## --- Tier 1 variants (strong clinical significance) ---

            ## A) Biomarker site _matches_ primary site of query tumor,
            ##    or pan-cancer biomarker (BM_PRIMARY_SITE = 'Any')
            ## B) Strong evidence - evidence levels A/B
            ## --> evidence item is tier-defining
            .data$ACTIONABILITY_TIER == 1 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              (.data$BM_PRIMARY_SITE == primary_site |
                 .data$BM_PRIMARY_SITE == "Any") &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$strong_regex) ~ "tier-defining",

            ## A) Biomarker site _matches_ primary site of query tumor
            ## B) Weak evidence - evidence levels C/D/E
            ## --> evidence item is providing additional support
            .data$ACTIONABILITY_TIER == 1 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              .data$BM_PRIMARY_SITE == primary_site &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex) ~ "additional",

            ## A) Biomarker site does not _match_ primary site of query tumor
            ##    and non pan-cancer biomarker (BM_PRIMARY_SITE != 'Any')
            ## B) Strong evidence - evidence levels A/B
            ## --> evidence item is providing additional support
            ##    (tier-1 was defined by a matching-site or pan-cancer strong item)
            .data$ACTIONABILITY_TIER == 1 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              (.data$BM_PRIMARY_SITE != primary_site &
                 .data$BM_PRIMARY_SITE != "Any") &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$strong_regex) ~ "additional",

            ## A) Biomarker site does not _match_ primary site of query tumor
            ##    and non pan-cancer biomarker (BM_PRIMARY_SITE != 'Any')
            ## B) Weak evidence - evidence levels C/D/E
            ## --> evidence item is providing additional support
            .data$ACTIONABILITY_TIER == 1 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              (.data$BM_PRIMARY_SITE != primary_site &
                 .data$BM_PRIMARY_SITE != "Any") &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex) ~ "additional",

            ## --- Tier 2 variants (potential clinical significance) ---

            ## A) Biomarker site does not _match_ primary site of query tumor
            ##    and non pan-cancer biomarker (BM_PRIMARY_SITE != 'Any')
            ## B) Strong evidence - evidence levels A/B
            ## --> evidence item is tier-defining
            .data$ACTIONABILITY_TIER == 2 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              (.data$BM_PRIMARY_SITE != primary_site &
                 .data$BM_PRIMARY_SITE != "Any") &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$strong_regex) ~ "tier-defining",

            ## A) Biomarker site _matches_ primary site of query tumor
            ## B) Weak evidence - evidence levels C/D/E
            ## --> evidence item is tier-defining
            .data$ACTIONABILITY_TIER == 2 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              .data$BM_PRIMARY_SITE == primary_site &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex) ~ "tier-defining",

            ## A) Biomarker site does not _match_ primary site of query tumor
            ##    and non pan-cancer biomarker (BM_PRIMARY_SITE != 'Any')
            ## B) Weak evidence - evidence levels C/D/E
            ## --> evidence item is providing additional support
            .data$ACTIONABILITY_TIER == 2 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              (.data$BM_PRIMARY_SITE != primary_site &
                 .data$BM_PRIMARY_SITE != "Any") &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL, bm_evidence$weak_regex) ~ "additional",

            ## --- Tier 3 variants (uncertain clinical significance) ---

            ## All weak evidence on tier-3 variants is tier-defining,
            ## regardless of pan-cancer or tumor-specific site — consistent
            ## with assign_bm_tier_support_ttagnostic, and correct now that
            ## pan-cancer weak contributes to tier-3 assignment
            .data$ACTIONABILITY_TIER == 3 &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering ~ "tier-defining",

            TRUE ~ as.character("other")
          )
        ) |>
        dplyr::mutate(BM_ACTIONABILITY_SUPPORT = factor(
          .data$BM_ACTIONABILITY_SUPPORT,
          levels = c("tier-defining",
                     "additional",
                     "additional-weak",
                     "other"),
          ordered = TRUE
        ))
    ) |>
      dplyr::arrange(
        .data$ACTIONABILITY_TIER,
        .data$BM_ACTIONABILITY_SUPPORT
      )
  }

  return(biomarker_items)

}

#' For all biomarker evidence items with specified confidence level,
#' assign whether these are tier-defining or providing additional
#' support - tumor-type agnostic query
#'
#' The function considers the strength of evidence (evidence levels) and
#' the match between biomarker site and primary site of query tumor.
#' Mark evidence considered tier-defining, while other biomarker evidence items
#' as additional support
#'
#' @param variants_tier_classified data frame with variant tier classifications
#' based on support from biomarker evidence item
#' @param vartype variant type (e.g. 'snv_indel', 'cna', 'fusion')
#' @param etype_for_tiering evidence type(s) to consider for tier classification
#' @param biomarker_items data frame with biomarker evidence items
#' @return data frame with biomarker evidence items classified as tier-defining
#' or providing additional support
#'
#' @export
#'
assign_bm_tier_support_ttagnostic <- function(
    variants_tier_classified = NULL,
    vartype = NULL,
    etype_for_tiering = c("predictive"),
    biomarker_items = NULL) {

  invisible(assertthat::assert_that(
    !is.null(variants_tier_classified) &
      is.data.frame(variants_tier_classified),
    msg = paste0("Argument var_df needs be of type data.frame")))
  invisible(
    assertthat::assert_that(
      !is.null(vartype) & is.character(vartype),
      vartype %in% c("snv_indel", "cna", "fusion"),
      msg = paste0("Argument 'vartype' needs to be one of
                   'snv_indel', 'cna' or 'fusion'"))
  )

  invisible(assertthat::assert_that(
    is.data.frame(biomarker_items),
    msg = paste0("Argument 'biomarker_items' needs be of type data.frame")))

  assertable::assert_colnames(
    variants_tier_classified,
    c("VAR_ID",
      "VARIANT_CLASS",
      "ACTIONABILITY",
      "ACTIONABILITY_TIER"),
    only_colnames = FALSE, quiet = TRUE)

  log4r_debug(paste0(
    "assign_bm_tier_support_ttagnostic - etype_for_tiering: ",
    paste(etype_for_tiering, collapse = ", ")))

  if (NROW(biomarker_items) > 0) {

    assertable::assert_colnames(
      biomarker_items,
      c("VAR_ID",
        "ENTREZGENE",
        "BM_EVIDENCE_LEVEL",
        "BM_PRIMARY_SITE"),
      only_colnames = FALSE, quiet = TRUE)

    ## Same fusion-specific join as in assign_bm_tier_support_ttspecific:
    ## use (VAR_ID, VARIANT_CLASS) only for fusions to avoid ENTREZGENE
    ## mismatches between single-gene DB entries and resolved pair values.
    tier_join_cols <- if (vartype == "fusion") {
      c("VAR_ID", "VARIANT_CLASS")
    } else {
      c("VAR_ID", "VARIANT_CLASS", "ENTREZGENE")
    }
    tier_lookup_cols <- c(tier_join_cols, "ACTIONABILITY_TIER")

    biomarker_items <- as.data.frame(
      biomarker_items |>
        dplyr::left_join(
          dplyr::select(
            variants_tier_classified,
            dplyr::all_of(tier_lookup_cols)
          ),
          by = tier_join_cols) |>
        dplyr::distinct() |>
        dplyr::select(
          c("VAR_ID",
            "VARIANT_CLASS",
            "ENTREZGENE",
            "ACTIONABILITY_TIER"),
          dplyr::everything()
        ) |>

        ## Now go through the max variant actionability tier associated with
        ## each evidence item and set them as tier-defining or supporting based
        ## on the strength of evidence and tumor type match
        dplyr::mutate(
          BM_ACTIONABILITY_SUPPORT = dplyr::case_when(

            ## A) Biomarker site is pan-cancer (tumor-agnostic)
            ## B) Strong evidence (CIViC A/B or OncoKB equivalents)
            ## C) Variant has a max tier of 1 (strong clinical significance)
            ## --> evidence item is tier-defining
            .data$BM_PRIMARY_SITE == "Any" &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              .data$ACTIONABILITY_TIER == 1 &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL,
                pcgrr::bm_evidence$strong_regex) ~ "tier-defining",

            ## A) Biomarker site is tumor-specific
            ## B) Strong evidence (CIViC A/B or OncoKB equivalents)
            ## C) Variant has a max tier of 2 (potential clinical significance)
            ## --> evidence item is tier-defining
            .data$BM_PRIMARY_SITE != "Any" &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              .data$ACTIONABILITY_TIER == 2 &
              stringr::str_detect(
                .data$BM_EVIDENCE_LEVEL,
                pcgrr::bm_evidence$strong_regex) ~ "tier-defining",

            ## A) Weak evidence (CIViC C/D/E or OncoKB equivalents)
            ## B) Variant has a max tier of 3 (uncertain clinical significance)
            ## --> evidence item is tier-defining
            stringr::str_detect(
              .data$BM_EVIDENCE_LEVEL,
              pcgrr::bm_evidence$weak_regex) &
              tolower(.data$BM_EVIDENCE_TYPE) %in% etype_for_tiering &
              .data$ACTIONABILITY_TIER == 3 ~ "tier-defining",

            !(tolower(.data$BM_EVIDENCE_TYPE) %in%
                etype_for_tiering) ~ "other",
            TRUE ~ as.character("additional")
          )
        ) |>
        dplyr::mutate(BM_ACTIONABILITY_SUPPORT = factor(
          .data$BM_ACTIONABILITY_SUPPORT,
          levels = c("tier-defining",
                     "additional",
                     "other"),
          ordered = TRUE
        ))
    ) |>
      dplyr::arrange(
        .data$ACTIONABILITY_TIER,
        .data$BM_ACTIONABILITY_SUPPORT
      )
  }

  return(biomarker_items)

}
