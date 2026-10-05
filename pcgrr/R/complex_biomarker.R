#' Properties of complex biomarkers (MSI-H, TMB-H)
#'
#' @param biomarker "MSI-H" or "TMB-H"
#' @return list with the variant type used in the biomarker data model
#' ('vartype'), the variant identifier ('var_id'), the name of the
#' corresponding molecular profile in CIViC ('profile_name'), the label used
#' for display ('label'), and the (pseudo) gene identifier used in the biomarker data model ('entrezgene')
#' @noRd
complex_biomarker_properties <- function(biomarker = "MSI-H") {
  props <- list(
    "MSI-H" = list(
      vartype = "msi",
      var_id = "MSI_1",
      profile_name = "MSI High",
      label = "MSI-High"),
    "TMB-H" = list(
      vartype = "tmb",
      var_id = "TMB_1",
      profile_name = "TMB High",
      label = "TMB-High"))
  if (!(length(biomarker) == 1 && biomarker %in% names(props))) {
    stop("Invalid biomarker - must be either 'MSI-H' or 'TMB-H'")
  }
  ## OncoKB annotates these biomarkers to the pseudo-gene 'Other Biomarkers'
  ## (Entrez gene identifier -2). The same identifier is used for evidence
  ## items from all sources, as they are not associated with a specific gene
  props[[biomarker]]$entrezgene <- "-2"
  props[[biomarker]]
}

#' Get biomarker match string for complex biomarker (MSI-H, TMB-H)
#'
#' Finds evidence items for the biomarker (single molecular profile named
#' "MSI High"/"TMB High" in CIViC/CGI, somatic origin), and formats them as
#' a biomarker match string, in the same format as used for biomarker matches
#' of other variant types (`<source>|<variant_id>|<evidence items>|<match type>`)
#'
#' @param biomarker "MSI-H" or "TMB-H"
#' @param ref_data PCGR reference data object
#' @return character (NA if there are no evidence items)
#' @noRd
get_complex_biomarker_match <- function(biomarker = "MSI-H", ref_data = NULL) {

  props <- complex_biomarker_properties(biomarker)
  clinical <- ref_data[["biomarker"]][["clinical"]]

  assertable::assert_colnames(
    clinical,
    c("BIOMARKER_SOURCE",
      "VARIANT_ID",
      "EVIDENCE_ID",
      "PRIMARY_SITE",
      "CANCER_TYPE",
      "CLINICAL_SIGNIFICANCE",
      "EVIDENCE_LEVEL",
      "EVIDENCE_TYPE",
      "VARIANT_ORIGIN",
      "MOLECULAR_PROFILE_NAME",
      "MOLECULAR_PROFILE_TYPE"),
    only_colnames = FALSE, quiet = TRUE)

  evidence <- clinical |>
    dplyr::filter(
      .data$MOLECULAR_PROFILE_NAME == props$profile_name &
        .data$MOLECULAR_PROFILE_TYPE == "Single" &
        !is.na(.data$VARIANT_ORIGIN) &
        .data$VARIANT_ORIGIN == "Somatic") |>
    dplyr::mutate(
      CLINICAL_SIGNIFICANCE_KEY = stringr::str_replace_all(
        stringr::str_replace_all(.data$CLINICAL_SIGNIFICANCE, " or ", "/"),
        "\\s+", "_"),
      EVIDENCE_LEVEL_KEY = stringr::str_replace(
        .data$EVIDENCE_LEVEL, ": .+", ""),
      EVIDENCE_TYPE_KEY = stringr::str_replace_all(
        .data$EVIDENCE_TYPE, "\\s+", "_"),
      PRIMARY_SITE_KEY = dplyr::case_when(
        !is.na(.data$PRIMARY_SITE) ~
          stringr::str_replace_all(.data$PRIMARY_SITE, "\\s+", "_"),
        .data$CANCER_TYPE %in%
          c("Cancer", "Solid Tumor", "Solid tumors") ~ "Any",
        TRUE ~ "Undefined"),
      EVIDENCE_KEY = paste(
        .data$EVIDENCE_ID, .data$PRIMARY_SITE_KEY,
        .data$CLINICAL_SIGNIFICANCE_KEY, .data$EVIDENCE_LEVEL_KEY,
        .data$EVIDENCE_TYPE_KEY, "Somatic", sep = ":"))

  if (NROW(evidence) == 0) {
    return(NA_character_)
  }

  match_type <- paste0("by_", props$vartype)
  matches <- evidence |>
    dplyr::group_by(.data$BIOMARKER_SOURCE, .data$VARIANT_ID) |>
    dplyr::summarise(
      EVIDENCE_ITEMS = paste(sort(unique(.data$EVIDENCE_KEY)), collapse = "&"),
      .groups = "drop") |>
    dplyr::mutate(
      MATCH = paste(
        .data$BIOMARKER_SOURCE, .data$VARIANT_ID,
        .data$EVIDENCE_ITEMS, match_type, sep = "|"))

  paste(matches$MATCH, collapse = ",")
}

#' Function that generates biomarker data for a complex biomarker
#' (MSI-H or TMB-H) for the PCGR report
#'
#' Evidence items are collected from CIViC/CGI (reference data) and from
#' OncoKB (web API, if enabled), and classified into tiers of clinical
#' significance (AMP/ASCO/CAP) for therapeutic sensitivity and resistance,
#' in the same manner as for other variant types (SNVs/InDels, CNAs, fusions).
#'
#' @param biomarker "MSI-H" or "TMB-H"
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#'
#' @return list with the biomarker record ('variant', 'variant_display') and
#' biomarker evidence ('bm_evidence'), i.e. a 'callset' for the biomarker
#'
#' @export
generate_report_data_complex_biomarker <- function(
    biomarker = "MSI-H",
    ref_data = NULL,
    settings = NULL) {

  props <- complex_biomarker_properties(biomarker)
  primary_site <- settings$conf$sample_properties$site
  if (is.null(primary_site) || is.na(primary_site)) {
    primary_site <- "Any"
  }

  log4r_info(paste0(
    "Mapping biomarker evidence - ", biomarker))

  callset <- list()
  callset[["variant"]] <- data.frame(
    SAMPLE_ID = settings$sample_id,
    VAR_ID = props$var_id,
    VARIANT_CLASS = props$vartype,
    ENTREZGENE = props$entrezgene,
    SAMPLE_ALTERATION = props$label,
    BIOMARKER_MATCH = get_complex_biomarker_match(biomarker, ref_data),
    stringsAsFactors = FALSE)
  callset[["variant_display"]] <- callset[["variant"]] |>
    dplyr::select(-c("BIOMARKER_MATCH"))
  callset[["bm_evidence"]] <- init_biomarker_content()

  oncokb_run <- isTRUE(as.numeric(settings$conf$oncokb$run) == 1)
  oncokb_exclusive <- isTRUE(as.numeric(settings$conf$oncokb$exclusive) == 1)

  ## CIViC/CGI - unless OncoKB is to be used exclusively
  if (!(oncokb_run && oncokb_exclusive) &&
      !is.na(callset[["variant"]]$BIOMARKER_MATCH)) {
    eitems <- map_biomarker_data(
      varcalls = callset[["variant"]],
      ref_data = ref_data,
      variant_origin = "Somatic",
      vartype = props$vartype)
    if (NROW(eitems) > 0) {
      callset[["bm_evidence"]][["eitems"]] <- eitems |>
        dplyr::mutate(ENTREZGENE = props$entrezgene)
    }
  }

  ## OncoKB
  if (oncokb_run) {
    okb_eitems <- fetch_oncokb_biomarker_eitems(
      biomarker = biomarker,
      oncotree_code = settings$conf$oncokb$oncotree_code,
      oncokb_token = settings$conf$oncokb$api_token)
    if (NROW(okb_eitems) > 0) {
      callset[["bm_evidence"]][["eitems"]] <- dplyr::bind_rows(
        callset[["bm_evidence"]][["eitems"]],
        okb_eitems)
    }
  }

  ## Tier classification (AMP/ASCO/CAP) - therapeutic sensitivity/resistance
  classified <- list()
  for (clnsig in c("therapeutic_sensitivity", "therapeutic_resistance")) {
    classified[[clnsig]] <- assign_amp_asco_cap_tiers(
      vartype = props$vartype,
      var_df = callset[["variant"]],
      primary_site = primary_site,
      clinical_significance = clnsig,
      biomarker_mapping_confidence = "medium",
      biomarker_items = callset[["bm_evidence"]][["eitems"]])

    for (elem in c("classification", "eitems")) {
      callset[["bm_evidence"]][[clnsig]][[elem]] <-
        classified[[clnsig]][["bm_evidence"]][[elem]]
    }
  }

  ## record with tier of therapeutic sensitivity
  if (NROW(classified[["therapeutic_sensitivity"]][["variant"]]) > 0) {
    callset[["variant"]] <-
      classified[["therapeutic_sensitivity"]][["variant"]]
    callset[["variant_display"]] <- callset[["variant"]] |>
      dplyr::select(-c("BIOMARKER_MATCH"))
    callset[["bm_evidence"]][["classification"]] <-
      classified[["therapeutic_sensitivity"]][["bm_evidence"]][["classification"]]
  }

  return(callset)
}
