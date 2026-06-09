#' Function that gathers data tables on actionable variants
#' for display in report (tier 1 + tier 2)
#'
#' @param rep report object
#' @param tier tier level(s) to consider for display (e.g. 1, 2)
#' @param etype_for_tiering evidence type(s) to consider for
#' tiering (e.g. 'predictive', 'prognostic', 'diagnostic')
#' @param clnsig clinical significance to consider for tiering
#' (e.g. 'therapeutic_sensitivity', 'therapeutic_resistance')
#' @param tier_defining_eitems_only consider only evidence items
#' that were used for tiering (e.g. for tier 1: only evidence items with
#' A-level evidence, for tier 2: only evidence items with B-level evidence).
#' If FALSE, all evidence items associated with each variant, not only
#' the tier-defining evidence items, will be considered for display in the report.
#' @param variant_category cna, snv_indel, or fusion
#'
#' @export
#'
prep_actble_display_tbl <- function(
    rep = NULL,
    tier = c(1,2),
    etype_for_tiering = c("predictive"),
    clnsig = "therapeutic_sensitivity",
    tier_defining_eitems_only = TRUE,
    variant_category = 'snv_indel') {

  if (is.null(rep)) {
    log4r_fatal("report object is NULL")
  }
  if (!is.numeric(tier) || !all(tier %in% c(1, 2, 3)) || length(tier) == 0) {
    log4r_fatal("tier must be a non-empty integer vector with values in {1, 2, 3}")
  }
  tier <- as.integer(unique(tier))
  if (!variant_category %in% names(rep$content)) {
    log4r_fatal(paste0(
      "rep$content object does not contain '", variant_category,"'"))
  }

  ## check variant_category is valid
  if (!variant_category %in% c("snv_indel", "cna", "fusion")) {
    log4r_fatal(
      "variant_category must be one of 'snv_indel', 'cna', or 'fusion'")
  }

  if (!"callset" %in% names(rep$content[[variant_category]])) {
    log4r_fatal("rep$content$variant_category object does not contain 'callset'")
  }

  callset <- rep$content[[variant_category]]$callset

  ## check that clinical_significance is one of ",therapeutic_sensitivity", "therapeutic_resistance",
  ## and make sure this is an element in the "bm_evidence" list object, which
  ## in turn should be a member of the "callset" list object for the variant class
  ## of interest (e.g. snv_indel, cna)
  if (!clnsig %in% c("therapeutic_sensitivity",
                     "therapeutic_resistance")) {
    log4r_fatal(
      paste0("clnsig must be one of ",
      "'therapeutic_resistance', 'therapeutic_sensitivity'"))
  }
  invisible(assertthat::assert_that(
    "bm_evidence" %in% names(callset),
    msg = paste0("rep$content$", variant_category,
                 "$callset object does not contain 'bm_evidence'"))
  )
   invisible(assertthat::assert_that(
     clnsig %in% names(callset$bm_evidence),
     msg = paste0("rep$content$", variant_category,
                  "$callset$bm_evidence object does not contain '",
                  clnsig,"'"))
   )

  eitems <-
    callset$bm_evidence[[clnsig]]$eitems
  var_classification <-
    callset$bm_evidence[[clnsig]]$classification

  ## For fusions the ENTREZGENE in biomarker evidence items (eitems) holds the
  ## single-gene DB value (e.g. "238" for ALK), whereas variant_display carries
  ## the resolved gene-pair (e.g. "238::238"). Joining on ENTREZGENE silently
  ## drops all fusion rows. Use (VAR_ID, VARIANT_CLASS) only for fusions;
  ## include ENTREZGENE for SNVs/InDels and CNAs as before.
  is_fusion_category <- isTRUE(variant_category == "fusion")
  display_join_cols <- if (is_fusion_category) {
    c("VAR_ID", "VARIANT_CLASS")
  } else {
    c("VAR_ID", "VARIANT_CLASS", "ENTREZGENE")
  }
  eitem_join_cols <- c("VAR_ID", "ACTIONABILITY_TIER", "VARIANT_CLASS",
                       if (!is_fusion_category) "ENTREZGENE")

  vars <- data.frame()

  if(clnsig == "therapeutic_sensitivity"){

    ## variants are by default organized for display in the
    ## report according to therapeutic sensitivity tiering.
    vars <-
      callset$variant_display |>
      dplyr::filter(
        !is.na(.data$ACTIONABILITY_TIER) &
          .data$ACTIONABILITY_TIER %in% tier)
  }else{
    ## For fusions, drop ENTREZGENE from variant_display before joining so
    ## that dplyr does not produce ENTREZGENE.x / ENTREZGENE.y suffixes
    ## (var_classification carries the canonical ENTREZGENE we want to keep).
    display_cols_to_drop <- c("ACTIONABILITY_TIER",
                              if (is_fusion_category) "ENTREZGENE")
    vars <-
      var_classification |>
      dplyr::filter(
        !is.na(.data$ACTIONABILITY_TIER) &
          .data$ACTIONABILITY_TIER %in% tier) |>
      dplyr::inner_join(
        dplyr::select(
          callset$variant_display,
          -dplyr::any_of(display_cols_to_drop)),
        by = display_join_cols
      )
  }

  if (NROW(vars) == 0) {
    log4r_info(
      paste0("No tier ", paste(tier, collapse = "/"), " variants found."))
    return(list(main = data.frame(), nested = data.frame()))
  }

  ## abundance column varies by variant category
  abundance_col <- switch(
    variant_category,
    snv_indel = "VAF_TUMOR",
    cna       = "CN_TOTAL",
    NULL
  )

  oncogenicity_col <- switch(
    variant_category,
    snv_indel = c("ONCOGENICITY",
                  "ONCOGENICITY_CODE"),
    NULL
  )

  hotspot_col <- switch(
    variant_category,
    snv_indel = "MUTATION_HOTSPOT",
    NULL
  )

  rctbl_recs <- list()
  rctbl_recs[['main']] <- data.frame()
  rctbl_recs[['nested']] <- data.frame()

  if (NROW(vars) > 0 & NROW(eitems) > 0) {
    ## debug: print colnames of vars and eitems

    log4r_debug(paste0(
      "Variant category: ", variant_category,
      " - colnames of vars: ", paste(colnames(vars), collapse = ", ")))

    ## When joining without ENTREZGENE (fusions), both sides carry an
    ## ENTREZGENE column. Drop it from vars so dplyr does not produce
    ## ENTREZGENE.x / ENTREZGENE.y suffixed columns; the single-gene
    ## ENTREZGENE from eitems is the one we want to keep downstream.
    vars_for_join <- if (is_fusion_category) {
      dplyr::select(vars, -dplyr::any_of("ENTREZGENE"))
    } else {
      vars
    }

    biomarker_var_eitems <- eitems |>
      dplyr::inner_join(
        vars_for_join,
        by = eitem_join_cols
      ) |>
      dplyr::distinct() |>
      dplyr::select(
        c("VAR_ID",
          "SAMPLE_ALTERATION",
          dplyr::any_of(oncogenicity_col),
          dplyr::any_of(abundance_col),
          dplyr::any_of(hotspot_col),
          "VARIANT_CLASS",
          "ENTREZGENE",
          "ACTIONABILITY_TIER",
          #"BM_EVIDENCE_ID",
          "BM_SOURCE_DB",
          "BM_REFERENCE",
          "BM_RATING",
          "BM_MOLECULAR_PROFILE",
          "BM_CANCER_TYPE",
          "BM_EVIDENCE_DESCRIPTION",
          "BM_EVIDENCE_TYPE",
          "BM_EVIDENCE_LEVEL",
          "BM_THERAPEUTIC_CONTEXT",
          "BM_CLINICAL_SIGNIFICANCE",
          "BM_PRIMARY_SITE",
          "BM_MAPPING_CONFIDENCE",
          "BM_RESOLUTION",
          "BM_ACTIONABILITY_SUPPORT")
      )

    ## only consider the evidence type that was used for tiering (e.g. predictive only,
    ## or all (predictive, prognostic, diagnostic)) also for display in the report
    if ("BM_EVIDENCE_TYPE" %in% colnames(biomarker_var_eitems) &
       NROW(biomarker_var_eitems) > 0) {
      biomarker_var_eitems <- biomarker_var_eitems |>
        dplyr::filter(
          tolower(.data$BM_EVIDENCE_TYPE) %in%
            etype_for_tiering)
    }

    if ("BM_ACTIONABILITY_SUPPORT" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0 &
        tier_defining_eitems_only == TRUE) {
      biomarker_var_eitems <-
        biomarker_var_eitems |>
        dplyr::filter(
          .data$BM_ACTIONABILITY_SUPPORT == "tier-defining")
    }

    if ("ONCOGENICITY" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0) {
      biomarker_var_eitems <- biomarker_var_eitems |>
        dplyr::mutate(
          ONCOGENICITY = stringr::str_replace_all(
            .data$ONCOGENICITY, "_", " "))
    }

    if (NROW(biomarker_var_eitems) > 0) {

      ## get the highest confidence level and resolution for
      ## each variant, based on all evidence items associated
      ## with the variant
      biomarker_top_resolution <-
        biomarker_var_eitems |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "BM_MAPPING_CONFIDENCE")
        ) |>
        dplyr::group_by(
          .data$ENTREZGENE,
          .data$VAR_ID,
        ) |>
        dplyr::summarise(
          BM_MAPPING_CONFIDENCE = paste(
            unique(sort(.data$BM_MAPPING_CONFIDENCE)),
            collapse = ","),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          BM_TOP_MAPPING_CONFIDENCE = dplyr::case_when(
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "high") ~ "high",
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "medium") ~ "medium",
            TRUE ~ "low"
          )
        ) |>
        dplyr::select(
          -c("BM_MAPPING_CONFIDENCE")) |>
        dplyr::distinct()

      ## across all evidence items, get the unique sources
      ## supporting actionability for each variant
      biomarker_source_support <-
        biomarker_var_eitems |>
        dplyr::group_by(
          .data$VAR_ID, .data$ENTREZGENE,
          .data$ACTIONABILITY_TIER
        ) |>
        dplyr::summarise(
          BM_SOURCES = paste(
            unique(sort(.data$BM_SOURCE_DB)),
            collapse = "|"),
          .groups = "drop") |>
        dplyr::distinct()

      ## for the main report table, we want to aggregate evidence items
      ## for each variant, and show the unique therapeutic contexts
      ## and primary sites associated with the variant, as well as the
      ## highest mapping confidence and resolution across all evidence items.
      ## For the nested table, we want to show all evidence items for each variant.
      ##
      rctbl_recs[['main']] <-
        biomarker_var_eitems |>
        dplyr::left_join(
          biomarker_top_resolution,
          by = c("VAR_ID","ENTREZGENE")
        ) |>
        dplyr::left_join(
          biomarker_source_support,
          by = c("VAR_ID","ENTREZGENE",
                 "ACTIONABILITY_TIER")
        ) |>
        dplyr::arrange(
          .data$BM_TOP_MAPPING_CONFIDENCE,
          dplyr::desc(.data$BM_RATING)
        ) |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "ACTIONABILITY_TIER",
            "SAMPLE_ALTERATION",
            "BM_SOURCES",
            dplyr::any_of(oncogenicity_col),
            dplyr::any_of(abundance_col),
            dplyr::any_of(hotspot_col),
            "BM_TOP_MAPPING_CONFIDENCE",
            "BM_THERAPEUTIC_CONTEXT",
            "BM_PRIMARY_SITE")
        ) |>
        dplyr::mutate(
          BM_THERAPEUTIC_CONTEXT = stringr::str_replace_all(
            .data$BM_THERAPEUTIC_CONTEXT, ",", ", "
          )
        ) |>
        dplyr::distinct() |>
        dplyr::group_by(
          dplyr::across(-c("BM_PRIMARY_SITE"))
        ) |>
        dplyr::summarise(
          BM_PRIMARY_SITE = paste(
            sort(.data$BM_PRIMARY_SITE),
            collapse = ", "),
          .groups = "drop") |>
        dplyr::ungroup() |>
        dplyr::mutate(
          BM_PRIMARY_SITE = stringr::str_replace_all(
            .data$BM_PRIMARY_SITE,
            "Any",
            "Multiple tumor types"
          )) |>
        dplyr::mutate(THERAPY_MATCH = paste0(
          "<b>",
          .data$BM_THERAPEUTIC_CONTEXT, "</b> (",
          .data$BM_PRIMARY_SITE, ")")) |>

        dplyr::group_by(
          dplyr::across(dplyr::all_of(
            c("VAR_ID",
              "ENTREZGENE",
              "ACTIONABILITY_TIER",
              "SAMPLE_ALTERATION",
              "BM_SOURCES",
              oncogenicity_col,
              abundance_col,
              hotspot_col,
              "BM_TOP_MAPPING_CONFIDENCE")
          ))
        ) |>
        dplyr::summarise(
          THERAPY_MATCH = paste(
            .data$THERAPY_MATCH, collapse=" | "),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          THERAPY_MATCH = paste0(
            " - ", .data$THERAPY_MATCH)
        ) |>
        dplyr::arrange(
          .data$ACTIONABILITY_TIER,
          .data$BM_TOP_MAPPING_CONFIDENCE,
        )

      rctbl_recs[['nested']] <- biomarker_var_eitems |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "ACTIONABILITY_TIER",
            "BM_THERAPEUTIC_CONTEXT",
            "BM_REFERENCE",
            "BM_MOLECULAR_PROFILE",
            "BM_EVIDENCE_LEVEL",
            "BM_SOURCE_DB",
            "BM_CANCER_TYPE",
            "BM_EVIDENCE_DESCRIPTION")
        ) |>
        dplyr::arrange(
          .data$BM_EVIDENCE_LEVEL
        ) |>
        dplyr::mutate(
          BM_THERAPEUTIC_CONTEXT = stringr::str_replace_all(
            .data$BM_THERAPEUTIC_CONTEXT, ",", ", "
          )
        ) |>
        dplyr::distinct()
    }
  }

  return(rctbl_recs)

}


#' Function that gathers data tables on prognostic variants
#' for display in report
#'
#' @param rep report object
#' @param tier tier level(s) to consider for display (e.g. 1, 2)
#' @param etype_for_tiering evidence type(s) to consider for
#' tiering (e.g. 'prognostic')
#' @param clnsig clinical significance to consider for tiering
#' (e.g. 'prognostic_better', 'prognostic_poor')
#' @param tier_defining_eitems_only consider only evidence items
#' that were used for tiering (e.g. for tier 1: only evidence items with
#' A-level evidence, for tier 2: only evidence items with B-level evidence).
#' If FALSE, all evidence items associated with each variant, not only
#' the tier-defining evidence items, will be considered for display in the report.
#' @param variant_category cna, snv_indel, or fusion
#'
#' @export
#'
prep_progn_display_tbl <- function(
    rep = NULL,
    tier = c(1,2),
    etype_for_tiering = c("prognostic"),
    clnsig = "prognostic_better",
    tier_defining_eitems_only = TRUE,
    variant_category = 'snv_indel') {

  if (is.null(rep)) {
    log4r_fatal("report object is NULL")
  }
  if (!is.numeric(tier) || !all(tier %in% c(1, 2)) || length(tier) == 0) {
    log4r_fatal("tier must be a non-empty integer vector with values in {1, 2}")
  }
  tier <- as.integer(unique(tier))
  if (!variant_category %in% names(rep$content)) {
    log4r_fatal(paste0(
      "rep$content object does not contain '", variant_category,"'"))
  }

  ## check variant_category is valid
  if (!variant_category %in% c("snv_indel", "cna", "fusion")) {
    log4r_fatal(
      "variant_category must be one of 'snv_indel', 'cna', or 'fusion'")
  }

  if (!"callset" %in% names(rep$content[[variant_category]])) {
    log4r_fatal("rep$content$variant_category object does not contain 'callset'")
  }

  callset <- rep$content[[variant_category]]$callset

  ## check that clinical_significance is one of "prognostic_better", "prognostic_poor",
  ## and make sure this is an element in the "bm_evidence" list object, which
  ## in turn should be a member of the "callset" list object for the variant class
  ## of interest (e.g. snv_indel, cna)
  if (!clnsig %in% c("prognostic_better", "prognostic_poor")) {
    log4r_fatal(
      paste0("clnsig must be one of ",
             "'prognostic_better', 'prognostic_poor'"))
  }
  invisible(assertthat::assert_that(
    "bm_evidence" %in% names(callset),
    msg = paste0("rep$content$", variant_category,
                 "$callset object does not contain 'bm_evidence'"))
  )
  invisible(assertthat::assert_that(
    clnsig %in% names(callset$bm_evidence),
    msg = paste0("rep$content$", variant_category,
                 "$callset$bm_evidence object does not contain '",
                 clnsig,"'"))
  )

  eitems <-
    callset$bm_evidence[[clnsig]]$eitems
  var_classification <-
    callset$bm_evidence[[clnsig]]$classification

  vars <-
    var_classification |>
    dplyr::filter(
      !is.na(.data$ACTIONABILITY_TIER) &
        .data$ACTIONABILITY_TIER %in% tier) |>
    dplyr::inner_join(
      dplyr::select(
        callset$variant_display,
        -c("ACTIONABILITY_TIER")),
      by = c("VAR_ID",
             "VARIANT_CLASS",
             "ENTREZGENE")
    )


  if (NROW(vars) == 0) {
    log4r_info(
      paste0("No tier ", paste(tier, collapse = "/"), " variants found."))
    return(list(main = data.frame(), nested = data.frame()))
  }

  ## abundance column varies by variant category
  abundance_col <- switch(
    variant_category,
    snv_indel = "VAF_TUMOR",
    cna       = "CN_TOTAL",
    NULL
  )

  oncogenicity_col <- switch(
    variant_category,
    snv_indel = "ONCOGENICITY",
    NULL
  )

  hotspot_col <- switch(
    variant_category,
    snv_indel = "MUTATION_HOTSPOT",
    NULL
  )

  rctbl_recs <- list()
  rctbl_recs[['main']] <- data.frame()
  rctbl_recs[['nested']] <- data.frame()

  if (NROW(vars) > 0 & NROW(eitems) > 0) {
    biomarker_var_eitems <- eitems |>
      dplyr::inner_join(
        vars,
        by = c("VAR_ID",
               "ACTIONABILITY_TIER",
               "VARIANT_CLASS",
               "ENTREZGENE")
      ) |>
      dplyr::distinct() |>
      dplyr::select(
        c("VAR_ID",
          "SAMPLE_ALTERATION",
          dplyr::any_of(oncogenicity_col),
          dplyr::any_of(abundance_col),
          dplyr::any_of(hotspot_col),
          "VARIANT_CLASS",
          "ENTREZGENE",
          "ACTIONABILITY_TIER",
          #"BM_EVIDENCE_ID",
          "BM_SOURCE_DB",
          "BM_REFERENCE",
          "BM_RATING",
          "BM_MOLECULAR_PROFILE",
          "BM_CANCER_TYPE",
          "BM_EVIDENCE_DESCRIPTION",
          "BM_EVIDENCE_TYPE",
          "BM_EVIDENCE_LEVEL",
          "BM_CLINICAL_SIGNIFICANCE",
          "BM_PRIMARY_SITE",
          "BM_MAPPING_CONFIDENCE",
          "BM_RESOLUTION",
          "BM_ACTIONABILITY_SUPPORT")
      )

    ## only consider the evidence type that was used for tiering (e.g. predictive only,
    ## or all (predictive, prognostic, diagnostic)) also for display in the report
    if ("BM_EVIDENCE_TYPE" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0) {
      biomarker_var_eitems <- biomarker_var_eitems |>
        dplyr::filter(
          tolower(.data$BM_EVIDENCE_TYPE) %in%
            etype_for_tiering)
    }

    if ("BM_ACTIONABILITY_SUPPORT" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0 &
        tier_defining_eitems_only == TRUE) {
      biomarker_var_eitems <-
        biomarker_var_eitems |>
        dplyr::filter(
          .data$BM_ACTIONABILITY_SUPPORT == "tier-defining")
    }

    if ("ONCOGENICITY" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0) {
      biomarker_var_eitems <- biomarker_var_eitems |>
        dplyr::mutate(
          ONCOGENICITY = stringr::str_replace_all(
            .data$ONCOGENICITY, "_", " "))
    }

    if (NROW(biomarker_var_eitems) > 0) {

      ## get the highest confidence level and resolution for
      ## each variant, based on all evidence items associated
      ## with the variant
      biomarker_top_resolution <-
        biomarker_var_eitems |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "BM_MAPPING_CONFIDENCE")
        ) |>
        dplyr::group_by(
          .data$ENTREZGENE,
          .data$VAR_ID,
        ) |>
        dplyr::summarise(
          BM_MAPPING_CONFIDENCE = paste(
            unique(sort(.data$BM_MAPPING_CONFIDENCE)),
            collapse = ","),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          BM_TOP_MAPPING_CONFIDENCE = dplyr::case_when(
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "high") ~ "high",
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "medium") ~ "medium",
            TRUE ~ "low"
          )
        ) |>
        dplyr::select(
          -c("BM_MAPPING_CONFIDENCE")) |>
        dplyr::distinct()


      ## across all evidence items, get the unique sources
      ## supporting actionability for each variant
      biomarker_source_support <-
        biomarker_var_eitems |>
        dplyr::group_by(
          .data$VAR_ID, .data$ENTREZGENE,
          .data$ACTIONABILITY_TIER
        ) |>
        dplyr::summarise(
          BM_SOURCES = paste(
            unique(sort(.data$BM_SOURCE_DB)),
            collapse = "|"),
          .groups = "drop") |>
        dplyr::distinct()

      ## for the main report table, we want to aggregate evidence items
      ## for each variant, and show the unique therapeutic contexts
      ## and primary sites associated with the variant, as well as the
      ## highest mapping confidence and resolution across all evidence items.
      ## For the nested table, we want to show all evidence items for each variant.
      ##
      rctbl_recs[['main']] <-
        biomarker_var_eitems |>
        dplyr::left_join(
          biomarker_top_resolution,
          by = c("VAR_ID","ENTREZGENE")
        ) |>
        dplyr::left_join(
          biomarker_source_support,
          by = c("VAR_ID","ENTREZGENE",
                 "ACTIONABILITY_TIER")
        ) |>
        dplyr::arrange(
          .data$BM_TOP_MAPPING_CONFIDENCE,
          dplyr::desc(.data$BM_RATING)
        ) |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "ACTIONABILITY_TIER",
            "SAMPLE_ALTERATION",
            "BM_SOURCES",
            dplyr::any_of(oncogenicity_col),
            dplyr::any_of(abundance_col),
            dplyr::any_of(hotspot_col),
            "BM_TOP_MAPPING_CONFIDENCE",
            "BM_CLINICAL_SIGNIFICANCE",
            "BM_PRIMARY_SITE")
        ) |>
        dplyr::distinct() |>
        dplyr::group_by(
          dplyr::across(-c("BM_PRIMARY_SITE"))
        ) |>
        dplyr::summarise(
          BM_PRIMARY_SITE = paste(
            sort(.data$BM_PRIMARY_SITE),
            collapse = ", "),
          .groups = "drop") |>
        dplyr::ungroup() |>
        dplyr::mutate(
          BM_PRIMARY_SITE = stringr::str_replace_all(
            .data$BM_PRIMARY_SITE,
            "Any",
            "Multiple tumor types"
          )) |>

        dplyr::mutate(PROGNOSTIC_OUTCOME = paste0(
          "<b>",
          .data$BM_CLINICAL_SIGNIFICANCE, "</b> (",
          .data$BM_PRIMARY_SITE, ")")) |>

        dplyr::group_by(
          dplyr::across(dplyr::all_of(
            c("VAR_ID",
              "ENTREZGENE",
              "ACTIONABILITY_TIER",
              "SAMPLE_ALTERATION",
              "BM_SOURCES",
              oncogenicity_col,
              abundance_col,
              hotspot_col,
              "BM_TOP_MAPPING_CONFIDENCE")
          ))
        ) |>
        dplyr::summarise(
          PROGNOSTIC_OUTCOME = paste(
            .data$PROGNOSTIC_OUTCOME, collapse=" | "),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          PROGNOSTIC_OUTCOME = paste0(
            " - ", .data$PROGNOSTIC_OUTCOME)
        ) |>
        dplyr::arrange(
          .data$ACTIONABILITY_TIER,
          .data$BM_TOP_MAPPING_CONFIDENCE,
        )

      rctbl_recs[['nested']] <- biomarker_var_eitems |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "BM_CLINICAL_SIGNIFICANCE",
            "ACTIONABILITY_TIER",
            "BM_REFERENCE",
            "BM_MOLECULAR_PROFILE",
            "BM_EVIDENCE_LEVEL",
            "BM_SOURCE_DB",
            "BM_CANCER_TYPE",
            "BM_EVIDENCE_DESCRIPTION")
        ) |>
        dplyr::arrange(
          .data$BM_EVIDENCE_LEVEL
        ) |>
        dplyr::distinct()
    }
  }

  return(rctbl_recs)

}


#' Function that gathers data tables on diagnostic variants
#'
#' @param rep report object
#' @param tier tier level(s) to consider for display (e.g. 1, 2)
#' @param etype_for_tiering evidence type(s) to consider for
#' tiering (e.g. 'diagnostic')
#' @param clnsig clinical significance to consider for tiering
#' (e.g. 'diagnostic_positive', 'diagnostic_negative')
#' @param tier_defining_eitems_only consider only evidence items
#' that were used for tiering (e.g. for tier 1: only evidence items with
#' A-level evidence, for tier 2: only evidence items with B-level evidence).
#' If FALSE, all evidence items associated with each variant, not only
#' the tier-defining evidence items, will be considered for display in the report.
#' @param variant_category cna, snv_indel, or fusion
#'
#' @export
#'
prep_diagn_display_tbl <- function(
    rep = NULL,
    tier = c(1,2),
    etype_for_tiering = c("diagnostic"),
    clnsig = "diagnostic_positive",
    tier_defining_eitems_only = TRUE,
    variant_category = 'snv_indel') {

  if (is.null(rep)) {
    log4r_fatal("report object is NULL")
  }
  if (!is.numeric(tier) || !all(tier %in% c(1, 2)) || length(tier) == 0) {
    log4r_fatal("tier must be a non-empty integer vector with values in {1, 2}")
  }
  tier <- as.integer(unique(tier))
  if (!variant_category %in% names(rep$content)) {
    log4r_fatal(paste0(
      "rep$content object does not contain '", variant_category,"'"))
  }

  ## check variant_category is valid
  if (!variant_category %in% c("snv_indel", "cna", "fusion")) {
    log4r_fatal(
      "variant_category must be one of 'snv_indel', 'cna', or 'fusion'")
  }

  if (!"callset" %in% names(rep$content[[variant_category]])) {
    log4r_fatal("rep$content$variant_category object does not contain 'callset'")
  }

  callset <- rep$content[[variant_category]]$callset

  ## check that clinical_significance is one of "prognostic_better", "prognostic_poor",
  ## and make sure this is an element in the "bm_evidence" list object, which
  ## in turn should be a member of the "callset" list object for the variant class
  ## of interest (e.g. snv_indel, cna)
  if (!clnsig %in% c("diagnostic_positive",
                     "diagnostic_negative")) {
    log4r_fatal(
      paste0("clnsig must be one of ",
             "'diagnostic_positive', 'diagnostic_negative'"))
  }
  invisible(assertthat::assert_that(
    "bm_evidence" %in% names(callset),
    msg = paste0("rep$content$", variant_category,
                 "$callset object does not contain 'bm_evidence'"))
  )
  invisible(assertthat::assert_that(
    clnsig %in% names(callset$bm_evidence),
    msg = paste0("rep$content$", variant_category,
                 "$callset$bm_evidence object does not contain '",
                 clnsig,"'"))
  )

  eitems <-
    callset$bm_evidence[[clnsig]]$eitems
  var_classification <-
    callset$bm_evidence[[clnsig]]$classification

  if (NROW(eitems) == 0 |
      NROW(var_classification) == 0) {
    log4r_info(
      paste0("No tier ", paste(tier, collapse = "/"), " variants found."))
    return(list(main = data.frame(), nested = data.frame()))
  }

  vars <-
    var_classification |>
    dplyr::filter(
      !is.na(.data$ACTIONABILITY_TIER) &
        .data$ACTIONABILITY_TIER %in% tier) |>
    dplyr::inner_join(
      dplyr::select(
        callset$variant_display,
        -c("ACTIONABILITY_TIER")),
      by = c("VAR_ID",
             "VARIANT_CLASS",
             "ENTREZGENE")
    )


  if (NROW(vars) == 0) {
    log4r_info(
      paste0("No tier ", paste(tier, collapse = "/"), " variants found."))
    return(list(main = data.frame(), nested = data.frame()))
  }

  ## abundance column varies by variant category
  abundance_col <- switch(
    variant_category,
    snv_indel = "VAF_TUMOR",
    cna       = "CN_TOTAL",
    NULL
  )

  oncogenicity_col <- switch(
    variant_category,
    snv_indel = "ONCOGENICITY",
    NULL
  )

  hotspot_col <- switch(
    variant_category,
    snv_indel = "MUTATION_HOTSPOT",
    NULL
  )

  rctbl_recs <- list()
  rctbl_recs[['main']] <- data.frame()
  rctbl_recs[['nested']] <- data.frame()

  if (NROW(vars) > 0 & NROW(eitems) > 0) {
    biomarker_var_eitems <- eitems |>
      dplyr::inner_join(
        vars,
        by = c("VAR_ID",
               "ACTIONABILITY_TIER",
               "VARIANT_CLASS",
               "ENTREZGENE")
      ) |>
      dplyr::distinct() |>
      dplyr::select(
        c("VAR_ID",
          "SAMPLE_ALTERATION",
          dplyr::any_of(oncogenicity_col),
          dplyr::any_of(abundance_col),
          dplyr::any_of(hotspot_col),
          "VARIANT_CLASS",
          "ENTREZGENE",
          "ACTIONABILITY_TIER",
          #"BM_EVIDENCE_ID",
          "BM_SOURCE_DB",
          "BM_REFERENCE",
          "BM_RATING",
          "BM_MOLECULAR_PROFILE",
          "BM_CANCER_TYPE",
          "BM_EVIDENCE_DESCRIPTION",
          "BM_EVIDENCE_TYPE",
          "BM_EVIDENCE_LEVEL",
          "BM_CLINICAL_SIGNIFICANCE",
          "BM_PRIMARY_SITE",
          "BM_MAPPING_CONFIDENCE",
          "BM_RESOLUTION",
          "BM_ACTIONABILITY_SUPPORT")
      )

    ## only consider the evidence type that was used for tiering (e.g. predictive only,
    ## or all (predictive, prognostic, diagnostic)) also for display in the report
    if ("BM_EVIDENCE_TYPE" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0) {
      biomarker_var_eitems <- biomarker_var_eitems |>
        dplyr::filter(
          tolower(.data$BM_EVIDENCE_TYPE) %in%
            etype_for_tiering)
    }

    if ("BM_ACTIONABILITY_SUPPORT" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0 &
        tier_defining_eitems_only == TRUE) {
      biomarker_var_eitems <-
        biomarker_var_eitems |>
        dplyr::filter(
          .data$BM_ACTIONABILITY_SUPPORT == "tier-defining")
    }

    if ("ONCOGENICITY" %in% colnames(biomarker_var_eitems) &
        NROW(biomarker_var_eitems) > 0) {
      biomarker_var_eitems <- biomarker_var_eitems |>
        dplyr::mutate(
          ONCOGENICITY = stringr::str_replace_all(
            .data$ONCOGENICITY, "_", " "))
    }

    if (NROW(biomarker_var_eitems) > 0) {

      ## get the highest confidence level and resolution for
      ## each variant, based on all evidence items associated
      ## with the variant
      biomarker_top_resolution <-
        biomarker_var_eitems |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "BM_MAPPING_CONFIDENCE")
        ) |>
        dplyr::group_by(
          .data$ENTREZGENE,
          .data$VAR_ID,
        ) |>
        dplyr::summarise(
          BM_MAPPING_CONFIDENCE = paste(
            unique(sort(.data$BM_MAPPING_CONFIDENCE)),
            collapse = ","),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          BM_TOP_MAPPING_CONFIDENCE = dplyr::case_when(
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "high") ~ "high",
            stringr::str_detect(
              .data$BM_MAPPING_CONFIDENCE, "medium") ~ "medium",
            TRUE ~ "low"
          )
        ) |>
        dplyr::select(
          -c("BM_MAPPING_CONFIDENCE")) |>
        dplyr::distinct()

      ## across all evidence items, get the unique sources
      ## supporting actionability for each variant
      biomarker_source_support <-
        biomarker_var_eitems |>
        dplyr::group_by(
          .data$VAR_ID, .data$ENTREZGENE,
          .data$ACTIONABILITY_TIER
        ) |>
        dplyr::summarise(
          BM_SOURCES = paste(
            unique(sort(.data$BM_SOURCE_DB)),
            collapse = "|"),
          .groups = "drop") |>
        dplyr::distinct()

      rctbl_recs[['main']] <-
        biomarker_var_eitems |>
        dplyr::left_join(
          biomarker_top_resolution,
          by = c("VAR_ID","ENTREZGENE")
        ) |>
        dplyr::left_join(
          biomarker_source_support,
          by = c("VAR_ID","ENTREZGENE",
                 "ACTIONABILITY_TIER")
        ) |>
        dplyr::arrange(
          .data$BM_TOP_MAPPING_CONFIDENCE,
          dplyr::desc(.data$BM_RATING)
        ) |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "ACTIONABILITY_TIER",
            "SAMPLE_ALTERATION",
            "BM_SOURCES",
            dplyr::any_of(oncogenicity_col),
            dplyr::any_of(abundance_col),
            dplyr::any_of(hotspot_col),
            "BM_TOP_MAPPING_CONFIDENCE",
            "BM_CLINICAL_SIGNIFICANCE",
            "BM_PRIMARY_SITE")
        ) |>
        dplyr::distinct() |>
        dplyr::group_by(
          dplyr::across(-c("BM_PRIMARY_SITE"))
        ) |>
        dplyr::summarise(
          BM_PRIMARY_SITE = paste(
            sort(.data$BM_PRIMARY_SITE),
            collapse = ", "),
          .groups = "drop") |>
        dplyr::ungroup() |>
        dplyr::mutate(
          BM_PRIMARY_SITE = stringr::str_replace_all(
            .data$BM_PRIMARY_SITE,
            "Any",
            "Multiple tumor types"
          )) |>

        dplyr::mutate(DIAGNOSTIC_EVIDENCE = paste0(
          "<b>",
          .data$BM_CLINICAL_SIGNIFICANCE, "</b> (",
          .data$BM_PRIMARY_SITE, ")")) |>

        dplyr::group_by(
          dplyr::across(dplyr::all_of(
            c("VAR_ID",
              "ENTREZGENE",
              "ACTIONABILITY_TIER",
              "SAMPLE_ALTERATION",
              "BM_SOURCES",
              oncogenicity_col,
              abundance_col,
              hotspot_col,
              "BM_TOP_MAPPING_CONFIDENCE")
          ))
        ) |>
        dplyr::summarise(
          DIAGNOSTIC_EVIDENCE = paste(
            .data$DIAGNOSTIC_EVIDENCE, collapse=" | "),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          DIAGNOSTIC_EVIDENCE = paste0(
            " - ", .data$DIAGNOSTIC_EVIDENCE)
        ) |>
        dplyr::arrange(
          .data$ACTIONABILITY_TIER,
          .data$BM_TOP_MAPPING_CONFIDENCE,
        )

      rctbl_recs[['nested']] <- biomarker_var_eitems |>
        dplyr::select(
          c("VAR_ID",
            "ENTREZGENE",
            "ACTIONABILITY_TIER",
            "BM_CLINICAL_SIGNIFICANCE",
            "BM_REFERENCE",
            "BM_MOLECULAR_PROFILE",
            "BM_EVIDENCE_LEVEL",
            "BM_SOURCE_DB",
            "BM_CANCER_TYPE",
            "BM_EVIDENCE_DESCRIPTION")
        ) |>
        dplyr::arrange(
          .data$BM_EVIDENCE_LEVEL
        ) |>
        dplyr::distinct()
    }
  }

  return(rctbl_recs)

}


#' Function that maps biomarker identifiers from VCF variant
#' annotation to full biomarker data tables for display in report
#'
#' @param varcalls variant calls data frame with biomarker identifiers
#' @param ref_data reference data object containing biomarker data
#' @param settings PCGR settings object
#' @param variant_origin variant origin for filtering biomarker evidence items
#' (e.g. "Somatic", "Germline", "Any")
#' @param vartype variant type for filtering biomarker evidence items
#' (e.g. "snv_indel", "cna", "fusion")
#' @return data frame with biomarker data for report display
#'
#' @export
#'
map_biomarker_data <- function(
    varcalls = NULL,
    ref_data = NULL,
    settings = NULL,
    variant_origin = "Somatic",
    vartype = "snv_indel") {


  assertthat::assert_that(
    !is.null(varcalls),
    msg = "Argument 'varcalls' must be provided")
  assertthat::assert_that(
    is.data.frame(varcalls),
    msg = "Argument 'varcalls' must be a valid data frame")

  ## check that ref_data object holds the required
  ## biomarker data tables

  assertthat::assert_that(
    !is.null(ref_data),
    msg = "Argument 'ref_data' must be provided")
  assertthat::assert_that(
    is.list(ref_data),
    msg = "Argument 'ref_data' must be a valid list object")
  assertthat::assert_that(
    "biomarker" %in% names(ref_data),
    msg = "Argument 'ref_data' must contain a 'biomarker' data table")
  assertthat::assert_that(
    "variant" %in% names(ref_data[['biomarker']]),
    msg = "Argument 'ref_data$biomarker' must contain a 'variant' data table")
  assertthat::assert_that(
    "clinical" %in% names(ref_data[['biomarker']]),
    msg = "Argument 'ref_data$biomarker' must contain a 'clinical' data table")
  assertthat::assert_that(
    "literature" %in% names(ref_data[['biomarker']]),
    msg = "Argument 'ref_data$biomarker' must contain a 'literature' data table")

  ## check for data frame structure of biomarker data tables
  assertthat::assert_that(
    is.data.frame(ref_data[['biomarker']][['variant']]),
    msg = "Argument 'ref_data$biomarker$variant' must be a valid data frame")

  assertthat::assert_that(
    is.data.frame(ref_data[['biomarker']][['clinical']]),
    msg = "Argument 'ref_data$biomarker$clinical' must be a valid data frame")

  assertthat::assert_that(
    is.data.frame(ref_data[['biomarker']][['literature']]),
    msg = "Argument 'ref_data$biomarker$literature' must be a valid data frame")

  ## check that literature, clinical and variant has the
  ## required columns for mapping, using assertable

  assertable::assert_colnames(
    ref_data[['biomarker']][['variant']],
    c("VARIANT_ID",
      "ENTREZGENE",
      "BIOMARKER_SOURCE"),
    only_colnames = FALSE, quiet = T
  )

  assertable::assert_colnames(
    ref_data[['biomarker']][['clinical']],
    c("EVIDENCE_ID",
      "CANCER_TYPE",
      "EVIDENCE_TYPE",
      "THERAPEUTIC_CONTEXT",
      "CLINICAL_SIGNIFICANCE",
      "DISEASE_ONTOLOGY_ID",
      "EVIDENCE_LEVEL",
      "EVIDENCE_DESCRIPTION",
      "MOLECULAR_PROFILE_NAME",
      "MOLECULAR_PROFILE_TYPE",
      "RATING",
      "EVIDENCE_DIRECTION"),
    only_colnames = FALSE, quiet = T
  )

  assertable::assert_colnames(
    ref_data[['biomarker']][['literature']],
    c("EVIDENCE_ID",
      "SOURCE_ID",
      "LINK",
      "NAME"),
    only_colnames = FALSE, quiet = T
  )

  biomarker_evidence_items <- data.frame()

  if ("BIOMARKER_MATCH" %in% colnames(varcalls) &
      "VAR_ID" %in% colnames(varcalls)) {

    ## find all sample variants that have been matched against
    ## biomarkers - (genomic/hgvsp/codon/exon/gene)
    biomarker_set <-
      varcalls |>
      dplyr::filter(!is.na(.data$BIOMARKER_MATCH))

    if (NROW(biomarker_set) > 0) {

      ## get citations of all biomarkers
      citations <- as.data.frame(
        ref_data[['biomarker']][['literature']] |>
          dplyr::select(
            c("EVIDENCE_ID",
              "SOURCE_ID",
              "LINK",
              "NAME")
          ) |>
          dplyr::mutate(
            CITATION = paste0(
              .data$SOURCE_ID,"|",
              .data$NAME
            )
          ) |>
          tidyr::separate_rows(
            c("EVIDENCE_ID"),
            sep=";"
          ) |>
          dplyr::group_by(
            .data$EVIDENCE_ID
          ) |>
          dplyr::summarise(
            CITATION = paste(
              unique(.data$CITATION), collapse = ", "
            ),
            CITATION_HTML = paste(
              unique(.data$LINK), collapse = ", "
            ),
            .groups = "drop"
          )
      )

      ## append more details of matched biomarkers from
      ## reference data on biomarkers (ref_dat$biomarker$variant)
      biomarker_evidence_items <-
        as.data.frame(
          biomarker_set |>
            dplyr::select(
              c("VAR_ID",
                "VARIANT_CLASS",
                "BIOMARKER_MATCH"),
            ) |>
            dplyr::distinct() |>
            tidyr::separate_rows(
              "BIOMARKER_MATCH", sep=",") |>
            tidyr::separate(
              "BIOMARKER_MATCH",
              into = c("BIOMARKER_SOURCE",
                       "VARIANT_ID",
                       "EVIDENCE_ITEMS",
                       "BIOMARKER_MATCHTYPE"),
              sep = "\\|",
              fill = "right"
            ) |>
            dplyr::mutate(
              VARIANT_ID = as.character(.data$VARIANT_ID)) |>
            dplyr::left_join(
              dplyr::select(
                ref_data[['biomarker']][['variant']],
                c("VARIANT_ID",
                  "ENTREZGENE",
                  "BIOMARKER_SOURCE")),
              by = c("VARIANT_ID",
                     "BIOMARKER_SOURCE"),
              relationship = "many-to-many") |>
            dplyr::rename("BIOMARKER_MATCH" = "BIOMARKER_MATCHTYPE") |>
            dplyr::mutate(BIOMARKER_RESOLUTION = dplyr::case_when(
              stringr::str_detect(
                .data$BIOMARKER_MATCH,"by_fusion") ~ "fusion",
              stringr::str_detect(
                .data$BIOMARKER_MATCH,"by_cna_segment") ~ "gene",
              stringr::str_detect(
                .data$BIOMARKER_MATCH,"by_genomic_coord") ~ "genomic",

              ## hgvsp match (no match at genomic level, but match
              ## at hgvsp level for principal transcript)
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_genomic_coord") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_hgvsp_principal") ~ "hgvsp",

              ## hgvsc match (no match at genomic level, but match at
              ## hgvsc level for principal transcript)
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_genomic_coord|by_hgvsp") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_hgvsc_principal") ~ "hgvsc",

              ## hgvsc match for non-principal transcript (no match at
              ## genomic level,
              ## but match at hgvsc level for non-principal transcript
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_genomic_coord|by_hgvsp") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_hgvsc_nonprincipal") ~ "hgvsc_nonprincipal",

              ## codon match (no match at genomic level, but match
              ## at codon level for principal
              ## transcript (no hgvsp match for principal transcript))
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_(genomic_coord|(hgvsp|hgvsc)_principal)") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_codon_principal") ~ "codon",

              ## Non-principal hgvsp match (no match at genomic level,
              ## but match at hgvsp level
              ## for non-principal transcript (no hgvsp match for
              ## principal transcript))
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_(genomic_coord|((hgvsp|hgvsc|codon)_principal))") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_hgvsp_nonprincipal") ~ "hgvsp_nonprincipal",

              ## Principal exon match (no match at genomic/hgvsp/codon level,
              ## but match at exon level
              ## for principal transcript)
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_(genomic_coord|((hgvsp|hgvsc|codon)_(principal|nonprincipal)))") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_exon_(mut|insertion|deletion)_principal") ~ "exon",

              ## Non-principal exon match (no match at genomic/hgvsp/codon level,
              ## but match at exon level
              ## for non-principal transcript (no exonmatch for principal transcript))
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_(genomic_coord|((hgvsp|hgvsc|codon)_(principal|nonprincipal)))") &
                !stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_exon_(mut|insertion|deletion)_principal") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "by_exon_(mut|insertion|deletion)_nonprincipal") ~ "exon_nonprincipal",

              ## gene region match (no match at genomic/hgvsp/codon/exon level,
              ## but match at aa region level)
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_(genomic_coord|hgvsp_|hgvsc_|codon_|exon_)") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,"by_aa_region") ~ "gene_region_mut",

              ## gene match (no match at genomic/hgvsp/codon/exon/aa region level,
              ##  but match at gene level)
              ## - here consider loss-of-function matches only,
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_(genomic_coord|hgvsp_|hgvsc_|codon_|exon_|aa_region_)") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,"by_gene_mut_lof") ~ "gene_lof_mut",

              ## gene match (no match at genomic/hgvsp/codon/exon/aa region level,
              ## but match at gene level)
              ## here consider non-LOF gene matches as well, since some non-LOF
              ## gene matches may be actionable (e.g. activating mutations in oncogenes)
              !stringr::str_detect(
                .data$BIOMARKER_MATCH,
                "by_(genomic_coord|hgvsp_|hgvsc_|codon_|exon_|aa_region_|gene_mut_lof)") &
                stringr::str_detect(
                  .data$BIOMARKER_MATCH,
                  "^by_gene_mut$") ~ "gene_mut",
              TRUE ~ as.character('other')
            )) |>

            tidyr::separate_rows(
              "EVIDENCE_ITEMS", sep="&") |>
            tidyr::separate(
              "EVIDENCE_ITEMS",
              into = c("EVIDENCE_ID",
                       "PRIMARY_SITE",
                       "CLINICAL_SIGNIFICANCE",
                       "EVIDENCE_LEVEL",
                       "EVIDENCE_TYPE",
                       "VARIANT_ORIGIN"),
              sep = ":"
            ) |>
            dplyr::select(
              -c("CLINICAL_SIGNIFICANCE",
                 "EVIDENCE_LEVEL",
                 "EVIDENCE_TYPE")) |>
            dplyr::mutate(PRIMARY_SITE = stringr::str_replace_all(
              .data$PRIMARY_SITE, "_"," "
            )) |>
            dplyr::left_join(
              citations, by = "EVIDENCE_ID",
              relationship = "many-to-many"
            ) |>
            dplyr::left_join(
              dplyr::select(
                ref_data[['biomarker']][['clinical']],
                c("EVIDENCE_ID",
                  "CANCER_TYPE",
                  "EVIDENCE_TYPE",
                  "THERAPEUTIC_CONTEXT",
                  "CLINICAL_SIGNIFICANCE",
                  "DISEASE_ONTOLOGY_ID",
                  "EVIDENCE_LEVEL",
                  "EVIDENCE_DESCRIPTION",
                  "MOLECULAR_PROFILE_NAME",
                  "MOLECULAR_PROFILE_TYPE",
                  "RATING",
                  "EVIDENCE_DIRECTION")
              ), by = c("EVIDENCE_ID"),
              relationship = "many-to-many"
            ) |>
            dplyr::rename(
              "BM_VARIANT_ID" = "VARIANT_ID",
              "BM_EVIDENCE_ID" = "EVIDENCE_ID",
              "BM_SOURCE_DB" = "BIOMARKER_SOURCE",
              "BM_RESOLUTION" = "BIOMARKER_RESOLUTION",
              "BM_MATCH" = "BIOMARKER_MATCH",
              "BM_PRIMARY_SITE" = "PRIMARY_SITE",
              "BM_EVIDENCE_TYPE" = "EVIDENCE_TYPE",
              "BM_CANCER_TYPE" = "CANCER_TYPE",
              "BM_DISEASE_ONTOLOGY_ID" = "DISEASE_ONTOLOGY_ID",
              "BM_VARIANT_ORIGIN" = "VARIANT_ORIGIN",
              "BM_EVIDENCE_LEVEL" = "EVIDENCE_LEVEL",
              "BM_EVIDENCE_DESCRIPTION" = "EVIDENCE_DESCRIPTION",
              "BM_THERAPEUTIC_CONTEXT" = "THERAPEUTIC_CONTEXT",
              "BM_CLINICAL_SIGNIFICANCE" = "CLINICAL_SIGNIFICANCE",
              "BM_CITATION" = "CITATION",
              "BM_REFERENCE" = "CITATION_HTML",
              "BM_RATING" = "RATING",
              "BM_EVIDENCE_DIRECTION" = "EVIDENCE_DIRECTION",
              "BM_MOLECULAR_PROFILE" = "MOLECULAR_PROFILE_NAME",
              "BM_MOLECULAR_PROFILE_TYPE" = "MOLECULAR_PROFILE_TYPE"
            ) |>
            dplyr::mutate(
              BM_RATING = dplyr::if_else(
                is.na(.data$BM_RATING),
                as.integer(2),
                as.integer(.data$BM_RATING)
              )
            ) |>
            dplyr::select(
              c("VAR_ID",
                "VARIANT_CLASS",
                "ENTREZGENE",
                "BM_SOURCE_DB",
                "BM_VARIANT_ID",
                "BM_EVIDENCE_ID",
                "BM_EVIDENCE_TYPE",
                "BM_EVIDENCE_LEVEL",
                "BM_EVIDENCE_DESCRIPTION",
                "BM_EVIDENCE_DIRECTION",
                "BM_CLINICAL_SIGNIFICANCE",
                "BM_VARIANT_ORIGIN",
                "BM_REFERENCE",
                "BM_CANCER_TYPE",
                "BM_DISEASE_ONTOLOGY_ID",
                "BM_PRIMARY_SITE",
                "BM_MATCH",
                "BM_RESOLUTION"),
              dplyr::everything()
            ) |>

            ## Make sure variant origin of patient/tumor
            ## matches reported variant origin of biomarker
            ## - Do not consider complex molecular profile types
            dplyr::filter(
              .data$BM_VARIANT_ORIGIN == variant_origin &
                .data$BM_MOLECULAR_PROFILE_TYPE == "Single") |>
            dplyr::distinct()
        )


      if (NROW(biomarker_evidence_items) > 0)
        biomarker_evidence_items <- biomarker_evidence_items |>
        dplyr::mutate(
          BM_MOLECULAR_PROFILE = dplyr::if_else(
            .data$BM_SOURCE_DB == "cgi",
            paste0(
              "<a href='https://www.cancergenomeinterpreter.org/2021/biomarkers'",
              " target='_blank'>",
              stringr::str_replace_all(
                .data$BM_MOLECULAR_PROFILE,
                ",",", "),
              "</a>"
            ),
            paste0(
              "<a href='https://civicdb.org/evidence/",
              stringr::str_replace(.data$BM_EVIDENCE_ID,"EID",""),
              "' target='_blank'>",
              stringr::str_replace_all(
                .data$BM_MOLECULAR_PROFILE,
                ",",", "),
              "</a>"
            )
          )
        ) |>
        dplyr::mutate(
          BM_MAPPING_CONFIDENCE = "low"
        ) |>
        dplyr::rename(
          BM_EVIDENCE_LEVEL_FULL = "BM_EVIDENCE_LEVEL"
        ) |>
        dplyr::mutate(
          BM_EVIDENCE_LEVEL = dplyr::case_when(
            startsWith(.data$BM_EVIDENCE_LEVEL_FULL,"A") ~ "A",
            startsWith(.data$BM_EVIDENCE_LEVEL_FULL,"B") ~ "B",
            startsWith(.data$BM_EVIDENCE_LEVEL_FULL,"C") ~ "C",
            startsWith(.data$BM_EVIDENCE_LEVEL_FULL,"D") ~ "D",
            startsWith(.data$BM_EVIDENCE_LEVEL_FULL,"E") ~ "E",
            TRUE ~ as.character("B")
          )
        )

      if (variant_origin == "Somatic" &
         NROW(biomarker_evidence_items) > 0) {

        if (vartype == "snv_indel") {
          ## For SNVs/indels with "somatic" variant origin, set a marker for
          ## "high" vs "low" confidence biomarker matches based on
          ## biomarker match type and variant oncogenicity/hotspot status, to
          ## prioritize display of biomarkers with higher confidence of association
          ## to the specific variant in the tumor, and to downweight display of
          ## biomarkers with lower confidence of association to the specific
          ## variant in the tumor, such as biomarkers matching at the "gene"
          ## level that are not supported by oncogenicity or hotspot evidence for
          ## the specific variant in the tumor. Do not filter, just set the mark

          ## by default (if oncogenicity/hotspot data not available),
          ## set "high" confidence for all, and then reset to "low" confidence
          ## for gene-level matches without genomic/hgvsp/codon/exon match
          ## evidence (i.e. gene-level matches without strong evidence of
          ## association to the specific variant in the tumor/patient)
          biomarker_evidence_items <- biomarker_evidence_items |>
            dplyr::mutate(
              BM_MAPPING_CONFIDENCE = dplyr::if_else(
                .data$BM_RESOLUTION == "hgvsp" |
                  .data$BM_RESOLUTION == "hgvsc" |
                  .data$BM_RESOLUTION == "genomic" |
                  .data$BM_RESOLUTION == "codon_principal" |
                  .data$BM_RESOLUTION == "exon_principal",
                "high",
                as.character(.data$BM_MAPPING_CONFIDENCE)
              )
            )

          if (NROW(varcalls) > 0 &
             "ONCOGENICITY" %in% colnames(varcalls) &
             "VAR_ID" %in% colnames(varcalls) &
             "ENTREZGENE" %in% colnames(varcalls) &
             "VARIANT_CLASS" %in% colnames(varcalls) &
             "MUTATION_HOTSPOT" %in% colnames(varcalls)) {

            biomarker_evidence_items$BM_MAPPING_CONFIDENCE <- NULL

            biomarker_evidence_items <-
              biomarker_evidence_items |>
              dplyr::left_join(
                dplyr::select(
                  varcalls,
                  c("VAR_ID",
                    "ONCOGENICITY",
                    "MUTATION_HOTSPOT",
                    "VARIANT_CLASS",
                    "ENTREZGENE")),
                by = c("VAR_ID",
                       "VARIANT_CLASS",
                       "ENTREZGENE")) |>
              dplyr::mutate(BM_MAPPING_CONFIDENCE = dplyr::case_when(
                .data$BM_RESOLUTION == "hgvsp" |
                  .data$BM_RESOLUTION == "genomic" |
                  .data$BM_RESOLUTION == "hgvsc" |
                  .data$BM_RESOLUTION == "codon" ~ "high",
                  .data$BM_RESOLUTION == "hgvsp_nonprincipal" |
                  .data$BM_RESOLUTION == "hgvsc_nonprincipal" |
                  .data$BM_RESOLUTION == "exon_nonprincipal" |
                  stringr::str_detect(
                    tolower(.data$ONCOGENICITY),"benign") ~ "low",
                (.data$BM_RESOLUTION == "gene_mut" |
                   .data$BM_RESOLUTION == "gene_lof_mut" |
                   .data$BM_RESOLUTION == "exon" |
                   .data$BM_RESOLUTION == "gene_region_mut") &
                  (stringr::str_detect(
                    tolower(.data$ONCOGENICITY),"oncogenic") |
                  !is.na(.data$MUTATION_HOTSPOT)) ~ "medium",
                TRUE ~ "low"
              )) |>
              dplyr::select(
                -c("ONCOGENICITY","MUTATION_HOTSPOT")
              ) |>
              dplyr::distinct()
          }
        }else{
          ## CNAs and fusions with "somatic" variant origin - set "high"
          ## confidence for all, since CNAs and fusions have
          ## gene-level resolution of biomarker matching
          biomarker_evidence_items <- biomarker_evidence_items |>
            dplyr::mutate(
              BM_MAPPING_CONFIDENCE = "high"
            )
        }
      }else{
        ## Germline variants
        if (NROW(biomarker_evidence_items) > 0 &
            variant_origin == "Germline" &
            "BM_RESOLUTION" %in% colnames(biomarker_evidence_items) &
            "BM_MAPPING_CONFIDENCE" %in% colnames(biomarker_evidence_items) &
            "CLASSIFICATION" %in% colnames(biomarker_evidence_items)){
          biomarker_evidence_items <- biomarker_evidence_items |>
            dplyr::mutate(
              BM_MAPPING_CONFIDENCE = dplyr::case_when(
                .data$BM_RESOLUTION == "hgvsp" |
                  .data$BM_RESOLUTION == "hgvsc" |
                  .data$BM_RESOLUTION == "genomic" |
                  .data$BM_RESOLUTION == "codon" ~ "high",
                ((.data$BM_RESOLUTION == "gene_mut" |
                   .data$BM_RESOLUTION == "gene_lof_mut" |
                   .data$BM_RESOLUTION == "exon" |
                   .data$BM_RESOLUTION == "gene_region_mut") &
                  (stringr::str_detect(
                    tolower(.data$CLASSIFICATION),
                    "pathogenic"))) ~ "medium",
                TRUE ~ "low"
              )
            )
        }
      }
    }

  }

  return(biomarker_evidence_items)

}



