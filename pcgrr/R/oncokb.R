#' Clean extracted evidence data from OncoKB
#'
#' @param evidence_df data.frame with OncoKB evidence
#' @param gene gene symbol
#' @param oncokb_root_gene root gene symbol from OncoKB annotation (for reference)
#' @param alteration name of alteration
#' @param vartype variant type (e.g., "snv_indel", "fusion", "cna")
#' @param profile_name molecular profile name (for display)
#' @param oncotree_code OncoTree code for tumor type
#'
#' @return Cleaned data.frame with standardized fields for
#' evidence level, clinical significance, and molecular profile formatting
#'
#' @export
#'
clean_oncokb_evidence <- function(
    evidence_df = NULL,
    gene = NA,
    oncokb_root_gene = NA,
    alteration = NA,
    vartype = NA,
    profile_name = NA,
    oncotree_code = NA){

  if (is.null(evidence_df) || NROW(evidence_df) == 0) {
    return(evidence_df)
  }

  gene_partners <- c()
  profile_name_url <- profile_name
  if(vartype == "fusion"){
    gene_partners <- unlist(
      stringr::str_split(gene, "-"))
    profile_name_url <-
      stringr::str_replace(
        profile_name, " - ", "%20")
  }

  assertable::assert_colnames(
    evidence_df,
    c("BM_EVIDENCE_LEVEL_FULL",
      "BM_RESOLUTION",
      "BM_REFERENCE",
      "BM_VARIANT_ORIGIN",
      "BM_CANCER_TYPE",
      "BM_MAPPING_CONFIDENCE",
      "BM_CLINICAL_SIGNIFICANCE",
      "BM_MOLECULAR_PROFILE"),
    only_colnames = FALSE,
    quiet = TRUE
  )

  evidence_df <- evidence_df |>
    dplyr::mutate(
      BM_MOLECULAR_PROFILE_ORIGINAL =
        BM_MOLECULAR_PROFILE
    ) |>
    tidyr::separate_rows(.data$BM_REFERENCE, sep=";") |>
    dplyr::mutate(
      BM_REFERENCE = dplyr::case_when(
        stringr::str_detect(.data$BM_REFERENCE, "^[0-9]{7,}$") ~
          paste0("<a href='https://pubmed.ncbi.nlm.nih.gov/",
                 .data$BM_REFERENCE, "' target='_blank'>",
                 .data$BM_REFERENCE, "</a>"),
        TRUE ~ as.character(.data$BM_REFERENCE)
      )
    ) |>
    dplyr::group_by(
      dplyr::across(-.data$BM_REFERENCE)) |>
      dplyr::summarise(
        BM_REFERENCE = paste(
          unique(.data$BM_REFERENCE), collapse = "; "),
        .groups = "drop") |>
    dplyr::mutate(BM_EVIDENCE_LEVEL = stringr::str_replace(
      .data$BM_EVIDENCE_LEVEL_FULL,"LEVEL_","")) |>
    dplyr::mutate(BM_CLINICAL_SIGNIFICANCE = dplyr::if_else(
      stringr::str_detect(
        .data$BM_EVIDENCE_LEVEL_FULL, "LEVEL_R"),
      "Resistance/Non-response",
      .data$BM_CLINICAL_SIGNIFICANCE
    )) |>
    dplyr::mutate(
      BM_MOLECULAR_PROFILE = dplyr::case_when(
        stringr::str_detect(
          .data$BM_MOLECULAR_PROFILE, "^Oncogenic") ~ paste0(
            gene, " - Oncogenic Mutations"),
        stringr::str_detect(
          .data$BM_MOLECULAR_PROFILE, "^[A-Z][0-9]{1,}( |$)") ~ paste0(
            gene, " - ", .data$BM_MOLECULAR_PROFILE),
        stringr::str_detect(
          .data$BM_MOLECULAR_PROFILE, "^Exon") ~ paste0(
            gene, " - ", stringr::str_replace_all(
              .data$BM_MOLECULAR_PROFILE,","," / ")),
        TRUE ~ profile_name
      )
    ) |>
    dplyr::mutate(
      BM_MOLECULAR_PROFILE = dplyr::case_when(
        vartype != "fusion" &
          .data$BM_VARIANT_ORIGIN == "Somatic" ~ paste0(
          "<a href='https://www.oncokb.org/gene/",gene,
          "/somatic/", alteration, "/",
          oncotree_code, "' target='_blank'>",
          .data$BM_MOLECULAR_PROFILE, "</a>"),
        vartype == "snv_indel" &
          .data$BM_VARIANT_ORIGIN == "Germline" ~ paste0(
            "<a href='https://www.oncokb.org/gene/",
            oncokb_root_gene,"/germline/Pathogenic%20Variants/",
            stringr::str_replace_all(
              stringr::str_replace_all(
              .data$BM_CANCER_TYPE, " ", "%20"), "/", "%2F"),
              "' target='_blank'>",
            .data$BM_MOLECULAR_PROFILE, "</a>"),
        vartype == "fusion" & length(gene_partners) == 2 ~ paste0(
          "<a href='https://www.oncokb.org/gene/",oncokb_root_gene,
          "/somatic/", profile_name_url, "/",
          oncotree_code, "' target='_blank'>",
          .data$BM_MOLECULAR_PROFILE, "</a>"),
        TRUE ~ as.character(.data$BM_MOLECULAR_PROFILE)
      )
    ) |>
    dplyr::mutate(BM_RESOLUTION = dplyr::case_when(
      stringr::str_detect(
        .data$BM_MOLECULAR_PROFILE, "Oncogenic Mutations|Exon ") &
        vartype == "snv_indel" ~ "gene_mut",
      TRUE ~ as.character(.data$BM_RESOLUTION)
    )) |>
    dplyr::mutate(BM_MAPPING_CONFIDENCE = dplyr::case_when(
      stringr::str_detect(
        .data$BM_MOLECULAR_PROFILE, "Oncogenic Mutations|Exon ") &
        vartype == "snv_indel" ~ "medium",
      TRUE ~ as.character(.data$BM_MAPPING_CONFIDENCE)
    ))

  return(evidence_df)

}



#' Extract therapeutic evidence items from OncoKB annotation
#'
#' @param oncokb_annotation OncoKB annotation list (from API)
#' @param gene Gene symbol (e.g., "BRAF")
#' @param alteration Alteration description
#' (e.g., HGVSp short format for SNVs/InDels)
#' @param vartype Variant type (e.g., "snv_indel", "fusion", "cna")
#' @param oncotree_code OncoTree code for tumor type (e.g., "BRCA" for breast cancer)
#' @param variant_id Variant identifier for tracking
#' @param match_by Matching strategy for SNVs/InDels (e.g., "hgvsp", "genomic")
#' @return Tibble with one row per treatment-cancer type combination
#'
#' @export
#'
extract_therapeutic_evidence <-
  function(oncokb_annotation,
           gene = NA,
           alteration = NA,
           vartype = NA,
           oncotree_code = NA,
           variant_id = NA,
           match_by = "hgvsp") {


    alteration2 <- stringr::str_replace(
      alteration,"p\\.","")
    profile_name <- paste(gene, alteration, sep = " - ")

    bmresolution = switch(
      vartype,
      "snv_indel" = match_by,
      "fusion" = "fusion",
      "cna" = "gene"
    )

    therapeutic_df <- data.frame()
    if (is.null(oncokb_annotation$treatments) ||
        length(oncokb_annotation$treatments) == 0) {
      return(therapeutic_df)
    }

    variant_origin <- "Somatic"
    if("pathogenic" %in% names(oncokb_annotation)){
      variant_origin <- "Germline"
    }

    # Helper function to extract first element from list/vector or return as-is
    extract_value <- function(x) {
      if (is.null(x)) return(NA)
      if (is.list(x) && length(x) == 1) return(x[[1]])
      if (is.vector(x) && length(x) == 1) return(x[1])
      if (is.vector(x) && length(x) > 1) return(paste(x, collapse = ","))
      if (is.list(x) && length(x) > 1) return(paste(sapply(x, as.character), collapse = ","))
      return(x)
    }

    # Extract each treatment as a row
    therapeutic_df <-
      purrr::map_df(oncokb_annotation$treatments, function(tx) {

        # Extract drug names
        drugs <- sapply(tx$drugs, function(d) d$drugName)
        # Sort the drug names alphabetically for consistent ordering
        drugs <- sort(drugs)

        drug_ncit_codes <- sapply(tx$drugs, function(d) {
          ifelse(is.null(d$ncitCode), NA, d$ncitCode)
        })

        # Extract cancer type - prioritize levelAssociatedCancerType
        cancer_type <- NA
        cancer_type_code <- NA
        cancer_type_main <- NA
        primary_site <- NA

        if (!is.null(tx$levelAssociatedCancerType)) {
          cancer_type <-
            extract_value(tx$levelAssociatedCancerType$name)
          cancer_type_code <-
            extract_value(tx$levelAssociatedCancerType$code)

          # If name is empty or NA, use mainType.name as fallback
          if (is.na(cancer_type) || cancer_type == "") {
            if (!is.null(tx$levelAssociatedCancerType$mainType)) {
              cancer_type <- extract_value(
                tx$levelAssociatedCancerType$mainType$name)
            }
          }

          # Store mainType separately for reference
          if (!is.null(tx$levelAssociatedCancerType$mainType)) {
            cancer_type_main <-
              extract_value(
                tx$levelAssociatedCancerType$mainType$name)
          }

          if (!is.null(tx$levelAssociatedCancerTyp$tissue)){
            primary_site <-
              extract_value(
                tx$levelAssociatedCancerType$tissue)
            if(primary_site == ""){
              primary_site <- "Any"
            }
            if(primary_site == "Bowel"){
              primary_site <- "Colon/Rectum"
            }
          }

        }

        # Create row for each cancer type - fit with existing biomarker
        # data model used in PCGR - one row per treatment-cancer type combination
        tibble::tibble(
          VAR_ID = variant_id,
          BM_SOURCE_DB = "oncokb",
          BM_VARIANT_ORIGIN = variant_origin,
          BM_RESOLUTION = bmresolution,
          BM_MAPPING_CONFIDENCE = "high",
          BM_EVIDENCE_TYPE = "Predictive",
          BM_EVIDENCE_DIRECTION = "Supports",
          BM_MOLECULAR_PROFILE = ifelse(
            is.null(tx$alterations) || length(tx$alterations) == 0,
            NA_character_,
            paste(tx$alterations, collapse = ", ")),
          BM_MOLECULAR_PROFILE_TYPE = "Single",
          BM_RATING = 5,
          BM_CLINICAL_SIGNIFICANCE = "Sensitivity/Response",
          BM_EVIDENCE_LEVEL_FULL = extract_value(tx$level),
          BM_FDA_LEVEL = extract_value(tx$fdaLevel),
          BM_THERAPEUTIC_CONTEXT = paste(drugs, collapse = ","),
          #BM_DRUG_CODE = paste(
            #drug_ncit_codes[!is.na(drug_ncit_codes)], collapse = ", "),
          BM_CANCER_TYPE = cancer_type,
          BM_PRIMARY_SITE = primary_site,
          BM_EVIDENCE_DESCRIPTION = extract_value(tx$description),
          BM_REFERENCE = if (!is.null(tx$pmids) && length(tx$pmids) > 0) {
            paste(sapply(tx$pmids, extract_value), collapse = ";")
          } else {
            NA
          },
        )
      })

    therapeutic_df <-
      clean_oncokb_evidence(
        evidence_df = therapeutic_df,
        oncokb_root_gene =
          oncokb_annotation$query$hugoSymbol,
        gene = gene,
        alteration = alteration2,
        vartype = vartype,
        oncotree_code = oncotree_code,
        profile_name = profile_name
      )

    return(therapeutic_df)
  }


#' Extract diagnostic implications from OncoKB annotation
#'
#' @param oncokb_annotation OncoKB annotation list
#' @param gene Gene symbol (e.g., "BRAF")
#' @param alteration Alteration description
#' (e.g., HGVSp short format for SNVs/InDels)
#' @param vartype Variant type (e.g., "snv_indel", "fusion", "cna")
#' @param oncotree_code OncoTree code for tumor type
#' (e.g., "BRCA" for breast cancer)
#' @param variant_id Variant identifier
#' @param match_by Matching strategy for SNVs/InDels (e.g., "hgvsp", "genomic")
#' @return Tibble with diagnostic evidence items
#' @export
extract_diagnostic_evidence <- function(
    oncokb_annotation = NULL,
    gene = NA,
    alteration = NA,
    vartype = NA,
    oncotree_code = NA,
    variant_id = NA,
    match_by = "hgvsp") {

  alteration2 <- stringr::str_replace(alteration, "p\\.", "")
  profile_name <- paste(gene, alteration, sep = " - ")

  diagnostic_df <- data.frame()
  if (is.null(oncokb_annotation$diagnosticImplications) ||
      length(oncokb_annotation$diagnosticImplications) == 0) {
    return(diagnostic_df)
  }

  variant_origin <- "Somatic"
  if("pathogenic" %in% names(oncokb_annotation)){
    variant_origin <- "Germline"
  }

  bmresolution = switch(
    vartype,
    "snv_indel" = match_by,
    "fusion" = "fusion",
    "cna" = "gene"
  )

  diagnostic_df <-
    purrr::map_df(oncokb_annotation$diagnosticImplications, function(dx) {
      tibble::tibble(
        VAR_ID = variant_id,
        BM_SOURCE_DB = "oncokb",
        BM_VARIANT_ORIGIN = variant_origin,
        BM_RESOLUTION = bmresolution,
        BM_MAPPING_CONFIDENCE = "high",
        BM_EVIDENCE_TYPE = "Diagnostic",
        BM_CLINICAL_SIGNIFICANCE = "Positive",
        BM_EVIDENCE_DIRECTION = "Supports",
        BM_MOLECULAR_PROFILE = ifelse(
          is.null(dx$alterations), NA, paste(dx$alterations, collapse = ", ")),
        BM_MOLECULAR_PROFILE_TYPE = "Single",
        BM_RATING = 5,
        BM_EVIDENCE_LEVEL_FULL = dx$level,
        BM_CANCER_TYPE = NA_character_,
        BM_EVIDENCE_DESCRIPTION = ifelse(
          is.null(dx$description), NA, dx$description),
        BM_REFERENCE = ifelse(
          is.null(dx$pmids), NA, paste(dx$pmids, collapse = ";"))
      )
    })

  diagnostic_df <-
    clean_oncokb_evidence(
      evidence_df = diagnostic_df,
      gene = gene,
      oncokb_root_gene =
        oncokb_annotation$query$hugoSymbol,
      alteration = alteration2,
      vartype = vartype,
      oncotree_code = oncotree_code,
      profile_name = profile_name
    )

  return(diagnostic_df)
}


#' Extract prognostic implications from OncoKB annotation
#'
#' @param oncokb_annotation OncoKB annotation list
#' @param gene Gene symbol (e.g., "BRAF")
#' @param alteration Alteration description
#' @param vartype Variant type (e.g., "snv_indel", "fusion", "cna")
#' (e.g., HGVSp short format for SNVs/InDels)
#' @param oncotree_code OncoTree code for tumor type
#' (e.g., "BRCA" for breast cancer)
#' @param variant_id Variant identifier
#' @param match_by Matching strategy for SNVs/InDels (e.g., "hgvsp", "genomic")
#' @return Tibble with prognostic evidence items
#' @export
extract_prognostic_evidence <- function(
    oncokb_annotation = NULL,
    gene = NA,
    alteration = NA,
    vartype = NA,
    oncotree_code = NA,
    variant_id = NA,
    match_by = "hgvsp") {

  alteration2 <- stringr::str_replace(alteration, "p\\.", "")
  profile_name <- paste(gene, alteration, sep = " - ")

  prognostic_df <- data.frame()
  if (is.null(oncokb_annotation$prognosticImplications) ||
      length(oncokb_annotation$prognosticImplications) == 0) {
    return(prognostic_df)
  }

  variant_origin <- "Somatic"
  if("pathogenic" %in% names(oncokb_annotation)){
    variant_origin <- "Germline"
  }

  bmresolution = switch(
    vartype,
    "snv_indel" = match_by,
    "fusion" = "fusion",
    "cna" = "gene"
  )

  prognostic_df <-
    purrr::map_df(oncokb_annotation$prognosticImplications, function(px) {
    tibble::tibble(
      VAR_ID = variant_id,
      BM_SOURCE_DB = "oncokb",
      BM_RESOLUTION = bmresolution,
      BM_MAPPING_CONFIDENCE = "high",
      BM_VARIANT_ORIGIN = variant_origin,
      BM_EVIDENCE_DIRECTION = "Supports",
      BM_EVIDENCE_TYPE = "Prognostic",
      BM_EVIDENCE_LEVEL_FULL = px$level,
      BM_CLINICAL_SIGNIFICANCE = "Poor Outcome",
      BM_MOLECULAR_PROFILE = ifelse(
        is.null(px$alterations), NA, paste(px$alterations, collapse = ", ")),
      BM_MOLECULAR_PROFILE_TYPE = "Single",
      BM_RATING = 5,
      BM_CANCER_TYPE = NA_character_,
      BM_EVIDENCE_DESCRIPTION = ifelse(
        is.null(px$description), NA, px$description),
      BM_REFERENCE = ifelse(
        is.null(px$pmids), NA, paste(px$pmids, collapse = ";")),
    )
  })

  prognostic_df <-
    clean_oncokb_evidence(
      evidence_df = prognostic_df,
      gene = gene,
      oncokb_root_gene =
        oncokb_annotation$query$hugoSymbol,
      alteration = alteration2,
      vartype = vartype,
      oncotree_code = oncotree_code,
      profile_name = profile_name
    )


  return(prognostic_df)
}


#' Extract mutation effect information
#'
#' param oncokb_annotation OncoKB annotation list
#' param variant_id Variant identifier
#' return Tibble with mutation effect data
#' export
# extract_mutation_effect <- function(oncokb_annotation, variant_id = NA) {
#
#   if (is.null(oncokb_annotation$mutationEffect)) {
#     return(NULL)
#   }
#
#   me <- oncokb_annotation$mutationEffect
#
#   # Extract citation PMIDs
#   pmids <- if (!is.null(me$citations) && length(me$citations) > 0) {
#     sapply(me$citations, function(c) {
#       if (!is.null(c$pmid)) c$pmid else NA
#     })
#   } else {
#     NA
#   }
#   pmids <- pmids[!is.na(pmids)]
#
#   tibble::tibble(
#     variant_id = variant_id,
#     mutation_effect = ifelse(
#       is.null(me$knownEffect), NA, me$knownEffect),
#     mutation_effect_description = ifelse(
#       is.null(me$description), NA, me$description),
#     mutation_effect_pmids = ifelse(
#       length(pmids) == 0, NA, paste(pmids, collapse = ";")),
#     n_mutation_effect_pmids = length(pmids)
#   )
# }


#' Extract oncogenicity classification
#'
#' param oncokb_annotation OncoKB annotation list
#' param variant_id Variant identifier
#' return Tibble with oncogenicity data
#' export
# extract_oncogenicity <- function(oncokb_annotation, variant_id = NA) {
#
#   tibble::tibble(
#     variant_id = variant_id,
#     oncogenic = ifelse(is.null(oncokb_annotation$oncogenic),
#                        NA, oncokb_annotation$oncogenic),
#     vus = ifelse(
#       is.null(oncokb_annotation$vus),
#       FALSE,
#       oncokb_annotation$vus),
#     highest_sensitive_level = ifelse(
#       is.null(oncokb_annotation$highestSensitiveLevel),
#       NA, oncokb_annotation$highestSensitiveLevel),
#     highest_resistance_level = ifelse(
#       is.null(oncokb_annotation$highestResistanceLevel),
#       NA, oncokb_annotation$highestResistanceLevel),
#     highest_diagnostic_level = ifelse(
#       is.null(oncokb_annotation$highestDiagnosticImplicationLevel),
#       NA, oncokb_annotation$highestDiagnosticImplicationLevel),
#     highest_prognostic_level = ifelse(
#       is.null(oncokb_annotation$highestPrognosticImplicationLevel),
#       NA, oncokb_annotation$highestPrognosticImplicationLevel)
#   )
# }


#' Extract clinical summaries
#'
#' param oncokb_annotation OncoKB annotation list
#' param variant_id Variant identifier
#' return Tibble with summary text fields
#' export
# extract_summaries <- function(oncokb_annotation, variant_id = NA) {
#
#   tibble::tibble(
#     variant_id = variant_id,
#     gene_summary = ifelse(
#       is.null(oncokb_annotation$geneSummary),
#       NA, oncokb_annotation$geneSummary),
#     variant_summary = ifelse(
#       is.null(oncokb_annotation$variantSummary),
#       NA, oncokb_annotation$variantSummary),
#     tumor_type_summary = ifelse(
#       is.null(oncokb_annotation$tumorTypeSummary),
#       NA, oncokb_annotation$tumorTypeSummary),
#     diagnostic_summary = ifelse(
#       is.null(oncokb_annotation$diagnosticSummary),
#       NA, oncokb_annotation$diagnosticSummary),
#     prognostic_summary = ifelse(
#       is.null(oncokb_annotation$prognosticSummary),
#       NA, oncokb_annotation$prognosticSummary)
#   )
# }


#' Extract complete structured annotation from OncoKB JSON
#'
#' @param oncokb_annotation OncoKB annotation list
#' @param gene gene symbol / gene fusion
#' @param alteration alteration description
#' (e.g., HGVSp short format for SNVs/InDels)
#' @param vartype variant type (e.g., "snv_indel", "fusion", "cna")
#' @param oncotree_code OncoTree code for tumor type
#' @param variant_id Variant identifier
#' @param match_by Matching strategy for SNVs/InDels (e.g., "hgvsp", "genomic")
#' @return Named list with all extracted components
#' @export
extract_complete_annotation <- function(
    oncokb_annotation,
    gene = NA,
    alteration = NA,
    vartype = NA,
    oncotree_code = NA,
    variant_id = NA,
    match_by = "hgvsp") {

  dplyr::bind_rows(
    extract_therapeutic_evidence(
      oncokb_annotation, gene, alteration,
      vartype, oncotree_code, variant_id, match_by),
    extract_diagnostic_evidence(
      oncokb_annotation, gene, alteration,
      vartype, oncotree_code, variant_id, match_by),
    extract_prognostic_evidence(
      oncokb_annotation, gene, alteration,
      vartype, oncotree_code, variant_id, match_by))

}


#' Fetch OncoKB annotation for SNV/InDel via protein change
#'
#' @param hugo_symbol Gene symbol (e.g., "BRAF")
#' @param protein_change Protein change in short format (e.g., "V600E")
#' @param oncotree_code OncoTree code/name (e.g., "SKIN", "BREAST")
#' @param oncokb_token OncoKB API token
#' @param base_api_url Optional base URL for OncoKB API (default: oncokb_base_api_url)
#' @param reference_genome Genome build, either "GRCh37" or "GRCh38" (default: "GRCh38")
#' @return List containing the complete JSON response from OncoKB API
#'
#' @export
#'
fetch_oncokb_hgvsp_annotation <-
  function(
    hugo_symbol = 'BRAF',
    protein_change = "V600E",
    oncotree_code = "SKIN",
    oncokb_token = NULL,
    base_api_url = NULL,
    reference_genome = "GRCh38") {

    # Validate inputs
    if (missing(oncokb_token) ||
        is.null(oncokb_token) ||
        oncokb_token == "") {
      stop("OncoKB token is required. Obtain from https://www.oncokb.org/account/settings")
    }

    # API endpoint
    base_url <- if (is.null(base_api_url)) {
      glue::glue("{oncokb_base_api_url}mutations/byProteinChange")
    } else {
      glue::glue("{base_api_url}mutations/byProteinChange")
    }

    # Build query parameters
    query_params <- list(
      hugoSymbol = hugo_symbol,
      alteration = protein_change,
      tumorType = oncotree_code,
      referenceGenome = reference_genome
    )

    if(is.null(oncotree_code) || length(oncotree_code) == 0 || is.na(oncotree_code)){
      query_params <- list(
        hugoSymbol = hugo_symbol,
        alteration = protein_change,
        referenceGenome = reference_genome
      )
    }

    # Make API request
    response <- tryCatch({
      httr::GET(
        url = base_url,
        query = query_params,
        httr::add_headers(
          Authorization = paste("Bearer", oncokb_token),
          Accept = "application/json"
        ),
        httr::timeout(30)
      )
    }, error = function(e) {
      warning(sprintf("API request failed for %s %s: %s",
                      hugo_symbol, protein_change, e$message))
      return(NULL)
    })

    # Check response status
    if (is.null(response)) {
      return(NULL)
    }

    if (httr::status_code(response) != 200) {
      warning(sprintf(
        "OncoKB API returned status %d for %s %s in %s: %s",
        httr::status_code(response),
        hugo_symbol,
        protein_change,
        oncotree_code,
        httr::content(response, "text", encoding = "UTF-8")
      ))
      return(NULL)
    }

    # Parse JSON response
    json_content <-
      httr::content(response, "text", encoding = "UTF-8")
    annotation <-
      jsonlite::fromJSON(json_content, simplifyVector = FALSE)

    return(annotation)
  }


#' Fetch OncoKB annotation for SNV/InDel via genomic change
#'
#' @param hgvsg Genomic change in HGVSg format (e.g., "7:g.140753336A>T")
#' @param oncotree_code Tumor type name
#' @param variant_origin somatic/germline
#' @param oncokb_token OncoKB API token
#' @param base_api_url Optional base URL for OncoKB API (default: oncokb_base_api_url)
#' @param reference_genome Genome build, either "GRCh37" or "GRCh38" (default: "GRCh38")
#' @return List containing the complete JSON response from OncoKB API
#'
#' @export
#'
fetch_oncokb_genomic_annotation <-
  function(
    hgvsg = "7:g.140753336A>T",
    oncotree_code = NULL,
    variant_origin = "somatic",
    oncokb_token = NULL,
    base_api_url = NULL,
    reference_genome = "GRCh38") {

  # Validate inputs
  if (missing(oncokb_token) ||
      is.null(oncokb_token) ||
      oncokb_token == "") {
    stop("OncoKB token is required. Obtain from https://www.oncokb.org/account/settings")
  }

  # check that variant_origin is either germline or somatic
  if (!variant_origin %in% c("somatic", "germline")) {
    stop("variant_origin must be either 'somatic' or 'germline'")
  }

  # API endpoint
  base_url <-
    glue::glue("{oncokb_base_api_url}mutations/byHGVSg")
  if(variant_origin == "germline"){
    base_url <-
      glue::glue("{oncokb_base_api_url}germline/mutations/byHGVSg")
  }
  if (!is.null(base_api_url)) {
    base_url <- glue::glue("{base_api_url}mutations/byHGVSg")
  }

  # Build query parameters
  query_params <- list(
    hgvsg = hgvsg,
    referenceGenome = reference_genome,
    tumorType = oncotree_code
  )

  if (is.null(oncotree_code) ||
      length(oncotree_code) == 0 ||
      is.na(oncotree_code)) {
    query_params <- list(
      hgvsg = hgvsg,
      referenceGenome = reference_genome
    )
  }

  # Make API request
  response <- tryCatch({
    httr::GET(
      url = base_url,
      query = query_params,
      httr::add_headers(
        Authorization = paste("Bearer", oncokb_token),
        Accept = "application/json"
      ),
      httr::timeout(30)
    )
  }, error = function(e) {
    warning(sprintf("API request failed for genomic change %s: %s",
                    hgvsg, e$message))
    return(NULL)
  })

  # Check response status
  if (is.null(response)) {
    return(NULL)
  }

  if (httr::status_code(response) != 200) {
    warning(sprintf(
      "OncoKB API returned status %d for genomic change %s in %s: %s",
      httr::status_code(response),
      hgvsg,
      oncotree_code,
      httr::content(response, "text", encoding = "UTF-8")
    ))
    return(NULL)
  }

  # Parse JSON response
  json_content <-
    httr::content(response, "text", encoding = "UTF-8")
  annotation <-
    jsonlite::fromJSON(json_content, simplifyVector = FALSE)

  return(annotation)
}


#' Fetch OncoKB annotation for gene fusion
#'
#' @param hugo_symbol_a First gene in fusion (5' partner)
#' @param hugo_symbol_b Second gene in fusion (3' partner)
#' @param oncotree_code OncoTree code/name (e.g., "SKIN", "BREAST")
#' @param oncokb_token OncoKB API token
#' @param base_api_url Optional base URL for OncoKB API (default: oncokb_base_api_url)
#' @param structural_variant_type Type of structural variant (default: "FUSION")
#' @param is_functional_fusion Whether fusion is functional (default: TRUE)
#' @return List containing the complete JSON response from OncoKB API
#' @export
fetch_oncokb_fusion_annotation <-
  function(
    hugo_symbol_a = "EML4",
    hugo_symbol_b = "ALK",
    oncotree_code = "LUNG",
    oncokb_token = NULL,
    base_api_url = NULL,
    structural_variant_type = "FUSION",
    is_functional_fusion = TRUE) {

    # Validate inputs
    if (missing(oncokb_token) ||
        is.null(oncokb_token) ||
        oncokb_token == "") {
      stop("OncoKB token is required")
    }

    # API endpoint
    base_url <- glue::glue(
      "{oncokb_base_api_url}structuralVariants")
    if (!is.null(base_api_url)) {
      base_url <- glue::glue("{base_api_url}structuralVariants")
    }


    # Build query parameters
    query_params <- list(
      hugoSymbolA = hugo_symbol_a,
      hugoSymbolB = hugo_symbol_b,
      structuralVariantType = structural_variant_type,
      isFunctionalFusion = tolower(as.character(is_functional_fusion)),
      tumorType = oncotree_code
    )

    if (is.null(oncotree_code) ||
          length(oncotree_code) == 0 ||
          is.na(oncotree_code)) {
      query_params <- list(
        hugoSymbolA = hugo_symbol_a,
        hugoSymbolB = hugo_symbol_b,
        structuralVariantType = structural_variant_type,
        isFunctionalFusion = tolower(as.character(is_functional_fusion))
      )
    }

    # Make API request
    response <- tryCatch({
      httr::GET(
        url = base_url,
        query = query_params,
        httr::add_headers(
          Authorization = paste("Bearer", oncokb_token),
          Accept = "application/json"
        ),
        httr::timeout(30)
      )
    }, error = function(e) {
      warning(sprintf("API request failed for %s-%s fusion: %s",
                      hugo_symbol_a, hugo_symbol_b, e$message))
      return(NULL)
    })

    # Check response
    if (is.null(response)) {
      return(NULL)
    }

    if (httr::status_code(response) != 200) {
      warning(sprintf(
        "OncoKB API returned status %d for %s-%s fusion in %s: %s",
        httr::status_code(response),
        hugo_symbol_a,
        hugo_symbol_b,
        oncotree_code,
        httr::content(response, "text", encoding = "UTF-8")
      ))
      return(NULL)
    }

    # Parse JSON response
    json_content <-
      httr::content(response, "text", encoding = "UTF-8")
    annotation <-
      jsonlite::fromJSON(json_content, simplifyVector = FALSE)

    return(annotation)
  }


#' Fetch OncoKB annotation for copy number alteration
#'
#' @param hugo_symbol Gene symbol
#' @param cna_type Type of CNA: "Amplification" or "Deletion"
#' @param oncotree_code OncoTree code/name (e.g., "SKIN", "BREAST")
#' @param oncokb_token OncoKB API token
#' @return List containing the complete JSON response from OncoKB API
#' @export
fetch_oncokb_cna_annotation <- function(
    hugo_symbol = NULL,
    cna_type = NULL,
    oncotree_code = NULL,
    oncokb_token = NULL) {

  # Validate inputs
  if (missing(oncokb_token) || is.null(oncokb_token) || oncokb_token == "") {
    log4r_fatal("OncoKB token is required")
  }

  # Map CNA type to OncoKB copyNameAlterationType
  cna_type_map <- c(
    "Amplification" = "AMPLIFICATION",
    "Deletion" = "DELETION"
  )

  if (!cna_type %in% names(cna_type_map)) {
    log4r_fatal(
      sprintf(
        paste0("Invalid CNA type: %s. Must be 'Amplification'",
               "or 'Deletion'", cna_type)))
  }

  # API endpoint
  base_url <- glue::glue(
    "{oncokb_base_api_url}copyNumberAlterations")

  # Build query parameters
  query_params <- list(
    hugoSymbol = hugo_symbol,
    copyNameAlterationType = cna_type_map[cna_type],
    tumorType = oncotree_code
  )

   if (is.null(oncotree_code) ||
       length(oncotree_code) == 0 ||
       is.na(oncotree_code)) {
     query_params <- list(
       hugoSymbol = hugo_symbol,
       copyNameAlterationType = cna_type_map[cna_type]
     )
   }

  # Make API request
  response <- tryCatch({
    httr::GET(
      url = base_url,
      query = query_params,
      httr::add_headers(
        Authorization = paste("Bearer", oncokb_token),
        Accept = "application/json"
      ),
      httr::timeout(30)
    )
  }, error = function(e) {
    warning(sprintf("API request failed for %s %s: %s",
                    hugo_symbol, cna_type, e$message))
    return(NULL)
  })

  # Check response
  if (is.null(response)) {
    return(NULL)
  }

  if (httr::status_code(response) != 200) {
    warning(sprintf(
      "OncoKB API returned status %d for %s %s in %s: %s",
      httr::status_code(response),
      hugo_symbol,
      cna_type,
      tumor_type,
      httr::content(response, "text", encoding = "UTF-8")
    ))
    return(NULL)
  }

  # Parse JSON response
  json_content <- httr::content(response, "text", encoding = "UTF-8")
  annotation <- jsonlite::fromJSON(json_content, simplifyVector = FALSE)

  return(annotation)
}


#' Process OncoKB MAF output files (both HGVSp and HGVSg) and fetch complete annotations
#'
#' @param maf_file_hgvsp Path to OncoKB-annotated MAF file using HGVSp (protein change) (optional)
#' @param maf_file_hgvsg Path to OncoKB-annotated MAF file using HGVSg (genomic change) (optional)
#' @param var_calls Data frame with variant calls (used for mapping of variant alteration)
#' @param oncotree_code OncoTree code used for OncoKB annotation
#' @param oncokb_token OncoKB API token
#' @param oncokb_base_api_url Optional base URL for OncoKB API (default: oncokb_base_api_url)
#' @param rate_limiting_delay Delay in seconds between API calls to respect rate limits
#' (default: 1 second)
#' @return Data frame with variant information and corresponding JSON annotations
#' @export
process_oncokb_maf <-
  function(
    maf_file_hgvsp = NULL,
    maf_file_hgvsg = NULL,
    var_calls = NULL,
    oncotree_code = NULL,
    oncokb_token = NULL,
    oncokb_base_api_url = NULL,
    rate_limiting_delay = 1) {

    # Validate that at least one MAF file is provided
    if (is.null(maf_file_hgvsp) && is.null(maf_file_hgvsg)) {
      log4r_fatal(
        "At least one MAF file (HGVSp or HGVSg) must be provided")
    }

    # Validate that var_calls is a non-null, a data.frame and non-empty
    # and containing 'VAR_ID' and 'ALTERATION' columns
    if (is.null(var_calls) || !is.data.frame(var_calls) ||
        nrow(var_calls) == 0 ||
        !all(c("VAR_ID", "ALTERATION","ENTREZGENE") %in%
             colnames(var_calls))) {
      log4r_fatal(
        paste0("var_calls must be a non-null data frame with ",
               "'VAR_ID','ENTREZGENE' and 'ALTERATION' columns"))
    }

    all_biomarker_evidence <- data.frame()
    all_variant_annotations <- data.frame()
    variant_keys <- character()  # Track unique variants to prevent duplicates

    # Process HGVSp file (protein changes)
    if (!is.null(maf_file_hgvsp) && file.exists(maf_file_hgvsp)) {
      maf_hgvsp <- readr::read_tsv(
        maf_file_hgvsp, show_col_types = FALSE)

      # Filter for variants in OncoKB genes regardless of ANNOTATED flag:
      # the MAF annotator batch endpoint misses variants covered only by gene-level
      # catchall entries (e.g. "Truncating Mutations"), while the individual API
      # endpoint correctly matches them — so we gate on GENE_IN_ONCOKB only and let
      # the API determine whether evidence exists
      maf_annotated_hgvsp <- maf_hgvsp |>
        dplyr::filter(
          .data$GENE_IN_ONCOKB == TRUE &
            .data$MUTATION_EFFECT != "Unknown")

      assertable::assert_colnames(
        maf_annotated_hgvsp,
        c("Hugo_Symbol",
          "Chromosome",
          "Tumor_Sample_Barcode",
          "Start_Position",
          "Reference_Allele",
          "Tumor_Seq_Allele2",
          "HGVSp_Short",
          "VAR_ID"),
        only_colnames = FALSE, quiet = TRUE
      )

      maf_annotated_hgvsp <- maf_annotated_hgvsp |>
        dplyr::rename(
          SYMBOL = "Hugo_Symbol",
          SAMPLE_ID = "Tumor_Sample_Barcode"
        ) |>
        dplyr::select(
          -dplyr::any_of(
            c("Chromosome", "Start_Position",
              "End_Position","Tumor_Seq_Allele1",
              "Reference_Allele", "Tumor_Seq_Allele2",
              "Variant_Classification",
              "HGVSp","NCBI_Build","ANNOTATED",
              "GENE_IN_ONCOKB","VARIANT_IN_ONCOKB")
          )
        ) |>
        dplyr::distinct() |>
        dplyr::filter(
          !is.na(HGVSp_Short) & HGVSp_Short != "") |>
        dplyr::select(
          SYMBOL, SAMPLE_ID, VAR_ID,
          dplyr::everything()
        ) |>
        dplyr::left_join(
          var_calls |>
            dplyr::select(
              VAR_ID,
              VARIANT_CLASS,
              ENTREZGENE,
              ALTERATION),
          by = "VAR_ID"
        )


      if (NROW(maf_annotated_hgvsp) > 0) {
        for (i in 1:NROW(maf_annotated_hgvsp)){
          var <- maf_annotated_hgvsp[i, ]

          if(var$VAR_ID %in% variant_keys) {
            next
          }

          all_variant_annotations <- dplyr::bind_rows(
            all_variant_annotations,
            var
          )

          ## log the symbol and hgvsp
          log4r_debug(
            sprintf("OncoKB web API retrieval - fetching annotation %d: %s %s",
                    length(variant_keys) + 1,
                    var$SYMBOL,
                    var$HGVSp_Short))

          # Fetch annotation via protein change
          annotation <- fetch_oncokb_hgvsp_annotation(
            hugo_symbol = var$SYMBOL,
            protein_change = var$HGVSp_Short,
            oncotree_code = oncotree_code,
            oncokb_token = oncokb_token,
            reference_genome = "GRCh38"
          )

          if (!is.null(annotation)) {

            evidence_df <- extract_complete_annotation(
              oncokb_annotation = annotation,
              gene = var$SYMBOL,
              oncotree_code = oncotree_code,
              alteration = var$ALTERATION,
              vartype = "snv_indel",
              variant_id = var$VAR_ID)

            if(NROW(evidence_df) > 0){

              evidence_df <- evidence_df |>
                dplyr::mutate(
                  VARIANT_CLASS = var$VARIANT_CLASS,
                  ENTREZGENE = var$ENTREZGENE,
                  BM_MATCH = "by_hgvsp_principal",
                ) |>
                dplyr::select(
                  VAR_ID,
                  VARIANT_CLASS,
                  ENTREZGENE,
                  BM_SOURCE_DB,
                  dplyr::everything()
                )

              all_biomarker_evidence <- dplyr::bind_rows(
                all_biomarker_evidence,
                evidence_df
              )
            }
            annotation <- NULL  # Clear annotation to save memory
            variant_keys <- c(variant_keys, var$VAR_ID)
          }
          # Rate limiting
          Sys.sleep(rate_limiting_delay)
        }
      }
    }

    # Process HGVSg file (considering non-protein changes)
    if (!is.null(maf_file_hgvsg) && file.exists(maf_file_hgvsg)) {
      maf_hgvsg <- readr::read_tsv(
        maf_file_hgvsg, show_col_types = FALSE)

      # Filter for variants in OncoKB genes without a protein change (HGVSg path),
      # regardless of ANNOTATED flag — same reasoning as HGVSp path above
      maf_annotated_hgvsg <- maf_hgvsg |>
        dplyr::filter(
          is.na(HGVSp_Short) | HGVSp_Short == "",
          .data$GENE_IN_ONCOKB == TRUE)

      assertable::assert_colnames(
        maf_annotated_hgvsg,
        c("Hugo_Symbol", "HGVSg",
          "Tumor_Sample_Barcode"),
        only_colnames = FALSE, quiet = TRUE
      )

      maf_annotated_hgvsg <- maf_annotated_hgvsg |>
        dplyr::rename(
          SYMBOL = Hugo_Symbol,
          SAMPLE_ID = Tumor_Sample_Barcode
        ) |>
        dplyr::select(
          -dplyr::any_of(
            c("Chromosome", "Start_Position",
              "End_Position","Tumor_Seq_Allele1",
              "Reference_Allele", "Tumor_Seq_Allele2",
              "Variant_Classification",
              "HGVSp","NCBI_Build","ANNOTATED",
              "GENE_IN_ONCOKB","VARIANT_IN_ONCOKB")
          )
        ) |>
        dplyr::distinct() |>
        dplyr::select(
          SYMBOL, SAMPLE_ID, VAR_ID,
          dplyr::everything()
        ) |>
        dplyr::left_join(
          var_calls |>
            dplyr::select(
              VAR_ID,
              VARIANT_CLASS,
              ENTREZGENE,
              ALTERATION),
          by = "VAR_ID"
        )


      if (NROW(maf_annotated_hgvsg) > 0) {
        for (i in 1:nrow(maf_annotated_hgvsg)) {
          var <- maf_annotated_hgvsg[i, ]
          if(var$VAR_ID %in% variant_keys) {
            next
          }

          all_variant_annotations <- dplyr::bind_rows(
            all_variant_annotations,
            var
          )

          ## log the symbol and hgvsp
          log4r_debug(
            sprintf("OncoKB web API retrieval - fetching annotation %d: %s %s",
                          length(variant_keys) + 1,
                          var$SYMBOL,
                          var$HGVSg))

          # Fetch annotation via genomic change
          annotation <- fetch_oncokb_genomic_annotation(
            hgvsg = var$HGVSg,
            oncotree_code = oncotree_code,
            variant_origin = "somatic",
            oncokb_token = oncokb_token,
            reference_genome = "GRCh38"
          )

          if (!is.null(annotation)) {

            evidence_df <- extract_complete_annotation(
              oncokb_annotation = annotation,
              gene = var$SYMBOL,
              oncotree_code = oncotree_code,
              alteration = var$ALTERATION,
              vartype = "snv_indel",
              variant_id = var$VAR_ID,
              match_by = "genomic")

            if(NROW(evidence_df) > 0){

              evidence_df <- evidence_df |>
                dplyr::mutate(
                  VARIANT_CLASS = var$VARIANT_CLASS,
                  ENTREZGENE = var$ENTREZGENE,
                  BM_MATCH = "by_genomic_coord"
                ) |>
                dplyr::select(
                  VAR_ID,
                  VARIANT_CLASS,
                  ENTREZGENE,
                  BM_SOURCE_DB,
                  dplyr::everything()
                )

              all_biomarker_evidence <- dplyr::bind_rows(
                all_biomarker_evidence,
                evidence_df
              )
            }
            annotation <- NULL  # Clear annotation to save memory
            variant_keys <- c(variant_keys, var$VAR_ID)
          }
          # Rate limiting
          Sys.sleep(rate_limiting_delay)
        }
      }
    }

    return(list('variants' = all_variant_annotations,
                'eitems' = all_biomarker_evidence))
  }


#' Process OncoKB fusion output file and fetch complete annotations
#'
#' @param fusion_file Path to OncoKB-annotated fusion file (TSV format)
#' @param var_calls Data frame with variant calls (used for mapping of variant alteration)
#' @param oncotree_code OncoTree code/name (e.g., "SKIN", "BREAST")
#' @param oncokb_token OncoKB API token
#' @param rate_limiting_delay Delay in seconds between API calls to respect rate limits
#' (default: 1 second)
#' @return List with fusion information and corresponding JSON annotations
#' @export
process_oncokb_fusion <-
  function(
    fusion_file = NULL,
    var_calls = NULL,
    oncotree_code = NULL,
    oncokb_token = NULL,
    rate_limiting_delay = 1) {


    all_biomarker_evidence <- data.frame()
    all_variant_annotations <- data.frame()

    assertthat::assert_that(
      !is.null(fusion_file) && file.exists(fusion_file),
      msg = "fusion_file must be a valid path to an existing file"
    )

    assertthat::assert_that(
      !is.null(var_calls) && is.data.frame(var_calls) &&
        all(c("VAR_ID","VARIANT_CLASS","ENTREZGENE") %in% colnames(var_calls)),
      msg = "var_calls must be a non-null data frame with 'VAR_ID' and 'ALTERATION' columns"
    )

    if (is.null(oncokb_token) || oncokb_token == "") {
      log4r_fatal(
        "OncoKB token is required. Obtain from https://www.oncokb.org/account/settings")
    }

    # Read OncoKB-annotated fusion file
    fusion_data <- readr::read_tsv(
      fusion_file, show_col_types = FALSE)

    assertable::assert_colnames(
      fusion_data,
      c("Fusion",
        "VAR_ID",
        "Tumor_Sample_Barcode",
        "ANNOTATED",
        "GENE_IN_ONCOKB",
        "VARIANT_IN_ONCOKB"),
      only_colnames = FALSE, quiet = TRUE
    )

    # Filter for annotated fusions
    annotated_fusions <- fusion_data |>
      dplyr::filter(.data$ANNOTATED == TRUE,
                    .data$GENE_IN_ONCOKB == TRUE,
                    .data$VARIANT_IN_ONCOKB == TRUE)

    if (nrow(annotated_fusions) == 0) {
      return(list('variants' = all_variant_annotations,
                  'eitems' = all_biomarker_evidence))
    }

    annotated_fusions <- annotated_fusions |>
      dplyr::rename(
        SAMPLE_ID = Tumor_Sample_Barcode) |>

      ## Join with var_calls to get VARIANT_CLASS and
      ## ENTREZGENE for each fusion
      dplyr::left_join(
        var_calls |>
          dplyr::select(
            VAR_ID,
            VARIANT_CLASS,
            ENTREZGENE),
        by = "VAR_ID"
      )

    results <- list()

    # Fetch annotations for each fusion
    for (i in 1:nrow(annotated_fusions)) {
      fusion <- annotated_fusions[i, ]

      # Parse fusion partners
      fusion_partners <- strsplit(fusion$Fusion, "-")[[1]]
      if (length(fusion_partners) != 2) {
        log4r_warn(
          sprintf("Cannot parse fusion: %s", fusion$Fusion))
        next
      }

      all_variant_annotations <- dplyr::bind_rows(
        all_variant_annotations,
        fusion
      )

      gene_a <- fusion_partners[1]
      gene_b <- fusion_partners[2]

      log4r_debug(
        sprintf("Fetching annotation %d/%d: %s-%s",
                i, nrow(annotated_fusions),
                gene_a, gene_b))

      # Fetch annotation
      annotation <- fetch_oncokb_fusion_annotation(
        hugo_symbol_a = gene_a,
        hugo_symbol_b = gene_b,
        oncotree_code = oncotree_code,
        oncokb_token = oncokb_token,
        is_functional_fusion = TRUE
      )

      if (!is.null(annotation)) {

        ## add each complex annotation object to the results object
        results[[i]] <- annotation

        evidence_df <- extract_complete_annotation(
          oncokb_annotation = annotation,
          gene = fusion$Fusion,
          oncotree_code = oncotree_code,
          alteration = "Fusion",
          vartype = "fusion",
          variant_id = fusion$VAR_ID)

        if(NROW(evidence_df) > 0){

          evidence_df <- evidence_df |>
            dplyr::mutate(
              VARIANT_CLASS = fusion$VARIANT_CLASS,
              ENTREZGENE = fusion$ENTREZGENE,
              BM_MATCH = "by_gene",
            ) |>
            dplyr::select(
              VAR_ID,
              VARIANT_CLASS,
              ENTREZGENE,
              BM_SOURCE_DB,
              dplyr::everything()
            )

          all_biomarker_evidence <- dplyr::bind_rows(
            all_biomarker_evidence,
            evidence_df
          )
        }
        annotation <- NULL  # Clear annotation to save memory
      }

      # Rate limiting
      Sys.sleep(rate_limiting_delay)
    }

    return(list('variants' = all_variant_annotations,
                'raw' = results,
                'eitems' = all_biomarker_evidence))
  }


#' Process OncoKB CNA output file and fetch complete annotations
#'
#' @param cna_file Path to OncoKB-annotated CNA file (TSV format)
#' @param oncotree_code OncoTree code/name (e.g., "SKIN", "BREAST")
#' @param oncokb_token OncoKB API token
#' @param var_calls Data frame with variant calls (used for mapping of variant alteration)
#' @param rate_limiting_delay Delay in seconds between API calls to respect rate
#' limits (default: 1 seconds)
#' @return List with CNA information and corresponding JSON annotations
#' @export
process_oncokb_cna <-
  function(
    cna_file,
    oncotree_code = NULL,
    oncokb_token = NULL,
    var_calls = NULL,
    rate_limiting_delay = 1) {


    assertthat::assert_that(
      !is.null(cna_file) && file.exists(cna_file),
      msg = "cna_file must be a valid path to an existing file"
    )

    assertthat::assert_that(
      !is.null(var_calls) && is.data.frame(var_calls) &&
        all(c("VAR_ID","VARIANT_CLASS","ENTREZGENE") %in% colnames(var_calls)),
      msg = "var_calls must be a non-null data frame with 'VAR_ID' and 'ALTERATION' columns"
    )

    if (is.null(oncokb_token) || oncokb_token == "") {
      log4r_fatal(
        "OncoKB token is required. Obtain from https://www.oncokb.org/account/settings")
    }

    # Read CNA file
    cna_data <- readr::read_tsv(
      cna_file, show_col_types = FALSE)

    all_biomarker_evidence <- data.frame()
    all_variant_annotations <- data.frame()

    assertable::assert_colnames(
      cna_data,
      c("Hugo_Symbol",
        "Tumor_Sample_Barcode",
        "VAR_ID",
        "Copy_Number_Alteration",
        "ANNOTATED",
        "GENE_IN_ONCOKB",
        "VARIANT_IN_ONCOKB"),
      only_colnames = FALSE, quiet = TRUE
    )

    # Filter for annotated CNAs
    annotated_cnas <- cna_data |>
      dplyr::filter(.data$ANNOTATED == TRUE,
                    .data$GENE_IN_ONCOKB == TRUE,
                    .data$VARIANT_IN_ONCOKB == TRUE)

    if (NROW(annotated_cnas) == 0) {
      return(list('variants' = all_variant_annotations,
                  'eitems' = all_biomarker_evidence))
    }

    # Fetch annotations for each CNA
    annotated_cnas <- annotated_cnas |>
      dplyr::rename(
        SAMPLE_ID = Tumor_Sample_Barcode) |>

      ## Join with var_calls to get VARIANT_CLASS and
      ## ENTREZGENE for each fusion
      dplyr::left_join(
        var_calls |>
          dplyr::select(
            VAR_ID,
            VARIANT_CLASS,
            SYMBOL,
            ENTREZGENE),
        by = c("VAR_ID" = "VAR_ID",
               "Hugo_Symbol" = "SYMBOL")
      )


    results <- list()

    for (i in 1:nrow(annotated_cnas)) {
      cna <- annotated_cnas[i, ]

      log4r_debug(
        sprintf("Fetching annotation %d/%d: %s %s",
                i, nrow(annotated_cnas),
                cna$Hugo_Symbol,
                cna$Copy_Number_Alteration))

      # Fetch annotation
      annotation <- fetch_oncokb_cna_annotation(
        hugo_symbol = cna$Hugo_Symbol,
        cna_type = cna$Copy_Number_Alteration,
        oncotree_code = oncotree_code,
        oncokb_token = oncokb_token
      )

      ## add each complex annotation object to the results object
      results[[i]] <- annotation

      all_variant_annotations <- dplyr::bind_rows(
        all_variant_annotations,
        cna
      )

      if (!is.null(annotation)) {

        evidence_df <- extract_complete_annotation(
          oncokb_annotation = annotation,
          gene = cna$Hugo_Symbol,
          oncotree_code = oncotree_code,
          vartype = "cna",
          alteration = stringr::str_to_title(
            cna$VARIANT_CLASS),
          variant_id = cna$VAR_ID)

        if(NROW(evidence_df) > 0){

          evidence_df <- evidence_df |>
            dplyr::mutate(
              VARIANT_CLASS = cna$VARIANT_CLASS,
              ENTREZGENE = cna$ENTREZGENE,
              BM_MATCH = "by_gene",
            ) |>
            dplyr::select(
              VAR_ID,
              VARIANT_CLASS,
              ENTREZGENE,
              BM_SOURCE_DB,
              dplyr::everything()
            )

          all_biomarker_evidence <- dplyr::bind_rows(
            all_biomarker_evidence,
            evidence_df
          )
        }
        annotation <- NULL  # Clear annotation to save memory
      }

      # Rate limiting
      Sys.sleep(rate_limiting_delay)
    }

    return(list('variants' = all_variant_annotations,
                'raw' = results,
                'eitems' = all_biomarker_evidence))
  }
