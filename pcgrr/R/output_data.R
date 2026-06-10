#' Build the SETTINGS sheet for the Excel workbook
#'
#' Returns a data frame with columns SECTION / PARAMETER / VALUE covering
#' the key PCGR run settings mirroring the HTML report Settings section.
#'
#' @param report PCGR report object
#' @return data.frame
#' @export
#'
get_settings_sheet <- function(report = NULL) {

  conf <- report$settings$conf
  s <- report$settings

  oncokb_active <-
    !is.null(conf$oncokb$api_token) &&
    !is.na(conf$oncokb$api_token) &&
    conf$oncokb$api_token != "None"
  bm_sources <- dplyr::case_when(
    oncokb_active && isTRUE(as.logical(conf$oncokb$exclusive)) ~ "OncoKB (exclusive)",
    oncokb_active ~ "CIViC + CGI + OncoKB",
    TRUE ~ "CIViC + CGI")

  rna_similarity <- if (isTRUE(as.logical(conf$expression$similarity_analysis))) "ON" else "OFF"

  rows <- list(
    ## General
    data.frame(SECTION = "General", PARAMETER = "PCGR version",
               VALUE = as.character(s$software$pcgr_version)),
    data.frame(SECTION = "General", PARAMETER = "Genome assembly",
               VALUE = as.character(s$genome_assembly)),
    data.frame(SECTION = "General", PARAMETER = "Biomarker sources",
               VALUE = bm_sources),
    data.frame(SECTION = "General", PARAMETER = "Actionability guidelines",
               VALUE = "AMP/ASCO/CAP"),
    data.frame(SECTION = "General", PARAMETER = "Oncogenicity guidelines",
               VALUE = "ClinGen/CGC/VICC")
  )

  ## SNV/InDel settings — only when SNV/InDel data present
  if (isTRUE(report$content$snv_indel$eval)) {
    as_conf  <- conf$somatic_snv$allelic_support
    sig_conf <- conf$somatic_snv$mutational_signatures
    rows <- c(rows, list(
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Show noncoding variants",
                 VALUE = as.character(conf$other$show_noncoding)),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Include germline findings",
                 VALUE = if (isTRUE(report$content$germline_classified$eval)) "ON" else "OFF"),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "MSI prediction",
                 VALUE = if (isTRUE(report$content$msi$eval)) "ON" else "OFF"),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Mutational burden estimation",
                 VALUE = if (isTRUE(report$content$tmb$eval)) "ON" else "OFF"),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "TMB algorithm",
                 VALUE = paste0("TMB_", conf$somatic_snv$tmb$tmb_display)),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Mutational signatures estimation",
                 VALUE = if (isTRUE(report$content$mutational_signatures$eval)) "ON" else "OFF"),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Signatures - min mutations required",
                 VALUE = as.character(sig_conf$mutation_limit)),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Signatures - all reference signatures",
                 VALUE = as.character(as.logical(sig_conf$all_reference_signatures))),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Signatures - include artefact signatures",
                 VALUE = as.character(as.logical(sig_conf$include_artefact_signatures))),
      data.frame(SECTION = "SNV/InDel", PARAMETER = "Signatures - min tumor-type prevalence (%)",
                 VALUE = as.character(sig_conf$prevalence_reference_signatures)),
      data.frame(SECTION = "VEP", PARAMETER = "Transcript set",
                 VALUE = if (isTRUE(as.logical(conf$vep$vep_gencode_basic)))
                   "GENCODE basic" else "GENCODE all"),
      data.frame(SECTION = "VEP", PARAMETER = "Transcript pick order",
                 VALUE = stringr::str_replace_all(conf$vep$vep_pick_order, ",", ", ")),
      data.frame(SECTION = "VEP", PARAMETER = "Regulatory regions annotation",
                 VALUE = as.character(as.logical(conf$vep$vep_regulatory))),
      data.frame(SECTION = "VEP", PARAMETER = "Buffer size",
                 VALUE = as.character(conf$vep$vep_buffer_size)),
      data.frame(SECTION = "VEP", PARAMETER = "Number of forks",
                 VALUE = as.character(conf$vep$vep_n_forks))
    ))
    ## Allelic support filters — only emit rows that are actually configured
    if (!is.null(as_conf$tumor_dp_tag) && as_conf$tumor_dp_tag != "_NA_" &&
        !is.null(as_conf$tumor_dp_min))
      rows <- c(rows, list(data.frame(SECTION = "SNV/InDel allelic filters",
        PARAMETER = "Min tumor depth (DP)", VALUE = as.character(as_conf$tumor_dp_min))))
    if (!is.null(as_conf$tumor_af_tag) && as_conf$tumor_af_tag != "_NA_" &&
        !is.null(as_conf$tumor_af_min))
      rows <- c(rows, list(data.frame(SECTION = "SNV/InDel allelic filters",
        PARAMETER = "Min tumor allelic fraction (AF)", VALUE = as.character(as_conf$tumor_af_min))))
    if (!is.null(as_conf$tumor_af_tag) && as_conf$tumor_af_tag != "_NA_" &&
        !is.null(as_conf$tumor_dp_tag) && as_conf$tumor_dp_tag != "_NA_" &&
        !is.null(as_conf$tumor_ad_min))
      rows <- c(rows, list(data.frame(SECTION = "SNV/InDel allelic filters",
        PARAMETER = "Min tumor allelic depth (AD)", VALUE = as.character(as_conf$tumor_ad_min))))
    if (!is.null(as_conf$control_dp_tag) && as_conf$control_dp_tag != "_NA_" &&
        !is.null(as_conf$control_dp_min))
      rows <- c(rows, list(data.frame(SECTION = "SNV/InDel allelic filters",
        PARAMETER = "Min control depth (DP)", VALUE = as.character(as_conf$control_dp_min))))
    if (!is.null(as_conf$control_af_tag) && as_conf$control_af_tag != "_NA_" &&
        !is.null(as_conf$control_af_max))
      rows <- c(rows, list(data.frame(SECTION = "SNV/InDel allelic filters",
        PARAMETER = "Max control allelic fraction (AF)", VALUE = as.character(as_conf$control_af_max))))
    if (!is.null(as_conf$control_dp_tag) && as_conf$control_dp_tag != "_NA_" &&
        !is.null(as_conf$control_af_tag) && as_conf$control_af_tag != "_NA_" &&
        !is.null(as_conf$control_ad_max))
      rows <- c(rows, list(data.frame(SECTION = "SNV/InDel allelic filters",
        PARAMETER = "Max control allelic depth (AD)", VALUE = as.character(as_conf$control_ad_max))))
  }

  ## CNA settings
  if (isTRUE(report$content$cna$eval)) {
    rows <- c(rows, list(
      data.frame(SECTION = "CNA", PARAMETER = "CNA thresholding mode",
                 VALUE = as.character(conf$somatic_cna$threshold_mode)),
      data.frame(SECTION = "CNA", PARAMETER = "Amplification threshold - absolute",
                 VALUE = as.character(conf$somatic_cna$amp_threshold_absolute)),
      data.frame(SECTION = "CNA", PARAMETER = "Amplification threshold - relative",
                 VALUE = as.character(conf$somatic_cna$amp_threshold_relative)),
      data.frame(SECTION = "CNA", PARAMETER = "Effective amplification threshold",
                 VALUE = as.character(conf$somatic_cna$amp_threshold_effective)),
      data.frame(SECTION = "CNA", PARAMETER = "Gain threshold - absolute",
                 VALUE = as.character(conf$somatic_cna$gain_threshold_absolute)),
      data.frame(SECTION = "CNA", PARAMETER = "Gain threshold - relative",
                 VALUE = as.character(conf$somatic_cna$gain_threshold_relative)),
      data.frame(SECTION = "CNA", PARAMETER = "Heterozygous deletion threshold - absolute",
                 VALUE = as.character(conf$somatic_cna$del_threshold_absolute)),
      data.frame(SECTION = "CNA", PARAMETER = "Heterozygous deletion threshold - relative",
                 VALUE = as.character(conf$somatic_cna$del_threshold_relative))
    ))
  }

  ## RNA fusion settings
  if (isTRUE(report$content$fusion$eval)) {
    rows <- c(rows, list(
      data.frame(SECTION = "RNA fusion", PARAMETER = "Min split reads required",
                 VALUE = as.character(conf$rna_fusion$min_split_reads))
    ))
  }

  ## Expression settings
  if (isTRUE(as.logical(report$content$expression$eval))) {
    rows <- c(rows, list(
      data.frame(SECTION = "Expression", PARAMETER = "RNA expression similarity analysis",
                 VALUE = rna_similarity)
    ))
  }

  dplyr::bind_rows(rows)
}

#' Build the DATA_VERSIONS sheet for the Excel workbook
#'
#' Returns a data frame with columns DATABASE / VERSION / DESCRIPTION / URL / LICENSE
#' covering all reference datasets used in the PCGR/CPSR run, mirroring the
#' "Dataset versions" callout in the HTML Documentation section.
#'
#' @param report PCGR or CPSR report object
#' @param wflow_pattern regex pattern to filter the wflow column
#'   (default `"pcgr"`; use `"cpsr"` for CPSR reports)
#' @param tool_label label shown in the bundle row DATABASE column
#'   (default `"PCGR"`)
#' @return data.frame
#' @export
#'
get_data_versions_sheet <- function(report = NULL,
                                    wflow_pattern = "pcgr",
                                    tool_label = "PCGR") {

  conf  <- report$settings$conf
  meta  <- report$settings$reference_data$source_metadata
  bundle_version <- as.character(report$settings$reference_data$version)

  oncokb_run <- isTRUE(as.integer(conf$oncokb$run) == 1L)

  rows <- list(
    data.frame(
      DATABASE    = paste0(tool_label, " data bundle"),
      VERSION     = bundle_version,
      DESCRIPTION = paste0("Reference data bundle for ", tool_label),
      URL         = "https://github.com/sigven/pcgr",
      LICENSE     = NA_character_,
      stringsAsFactors = FALSE)
  )

  seen <- character(0)

  for (i in seq_len(NROW(meta))) {
    wflow <- meta[i, "wflow"]
    if (!stringr::str_detect(wflow, wflow_pattern)) next

    source      <- meta[i, "source_abbreviation"]
    source_full <- meta[i, "source"]
    source_type <- meta[i, "source_type"]
    description <- meta[i, "source_description"]
    url         <- meta[i, "source_url"]
    version     <- meta[i, "source_version"]
    license     <- meta[i, "source_license"]

    ## Same deduplication rules as documentation.qmd
    if (source == "dbnsfp"         && source_type == "gene")      next
    if (source == "depmap"         && source_type != "other")     next
    if (source == "tcga_pancan_2018" && source_type == "phenotype") next

    ## Skip OncoKB row if not in use
    if (source == "oncokb" && !oncokb_run) next

    ## Each source only once
    if (source %in% seen) next
    seen <- c(seen, source)

    ## OncoKB: override version with the one stored in conf
    if (source == "oncokb" && oncokb_run) {
      version <- paste0(
        conf$oncokb$data_version, " - ", conf$oncokb$data_release_date)
    }

    rows <- c(rows, list(data.frame(
      DATABASE    = source_full,
      VERSION     = if (is.na(version)) NA_character_ else as.character(version),
      DESCRIPTION = as.character(description),
      URL         = as.character(url),
      LICENSE     = as.character(license),
      stringsAsFactors = FALSE)))
  }

  dplyr::bind_rows(rows)
}

#' Function that produces the contents of sheets for an Excel report
#' of PCGR output
#'
#' @param report PCGR report object
#'
#' @export
#'
get_excel_sheets <- function(report = NULL) {

  if (is.null(report)) {
    stop("report must be provided")
  }

  if (!((base::typeof(report$content) == "list") == TRUE)) {
    stop("report$content must be a list")
  }

  ## Settings and data versions sheets (first two in workbook)
  excel_sheets <- list()
  excel_sheets[['SETTINGS']] <- get_settings_sheet(report)
  excel_sheets[['DATA_VERSIONS']] <- get_data_versions_sheet(report)

  ## Sample and assay properties
  excel_sheets[['SAMPLE_ASSAY']] <- data.frame()
  sample_id <- report$settings$sample_id

  if (!((base::typeof(report$content$sample_properties) == "list") == TRUE)) {
    stop("report$content$sample_properties must be a list")
  }

  callsets <- list()
  callsets[['snv_indel']] <- report$content$snv_indel$callset
  callsets[['cna']] <- report$content$cna$callset
  callsets[['fusion']] <- report$content$fusion$callset

  .empty_alt <- data.frame(
    VAR_ID = character(),
    ENTREZGENE = integer(),
    VARIANT_CLASS = character(),
    SAMPLE_ALTERATION = character(),
    stringsAsFactors = FALSE)

  sample_alteration <- list(
    snv_indel = .empty_alt,
    cna       = .empty_alt,
    fusion    = .empty_alt)

  for(t in c('snv_indel','cna','fusion')) {
    if (!is.null(callsets[[t]])) {
      if (NROW(callsets[[t]]$variant) > 0) {
        sample_alteration[[t]] <-
          callsets[[t]]$variant |>
          dplyr::select(
            "VAR_ID",
            "ENTREZGENE",
            "VARIANT_CLASS",
            "SAMPLE_ALTERATION") |>
          dplyr::distinct()
      }
    }
  }

  ## initialize sample and assay information
  sample_assay <-
    dplyr::bind_rows(
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'SAMPLE',
        PROPERTY = 'SITE', VALUE = NA
      ),
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'SAMPLE',
        PROPERTY = 'SEX', VALUE = NA
      ),
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'SAMPLE',
        PROPERTY = 'TUMOR_PURITY', VALUE = NA
      ),
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'SAMPLE',
        PROPERTY = 'TUMOR_PLOIDY', VALUE = NA
      ),
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'ASSAY',
        PROPERTY = 'TYPE', VALUE = NA
      ),
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'ASSAY',
        PROPERTY = 'MODE', VALUE = NA
      ),
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'ASSAY',
        PROPERTY = 'EFFECTIVE_TARGET_SIZE_MB', VALUE = NA
      ),
    )

  for (elem in c('SITE','SEX','TUMOR_PURITY','TUMOR_PLOIDY')) {
    if (tolower(elem) %in% names(report$content$sample_properties)) {
        sample_assay[sample_assay$PROPERTY == elem, 'VALUE'] <-
          report$content$sample_properties[[tolower(elem)]]
      }
    }

  for (elem in c('TYPE','MODE','EFFECTIVE_TARGET_SIZE_MB')) {
    if (tolower(elem) %in% names(report$content$assay_properties)) {
      sample_assay[sample_assay$PROPERTY == elem, 'VALUE'] <-
        report$content$assay_properties[[tolower(elem)]]
    }
  }

  excel_sheets[['SAMPLE_ASSAY']] <- sample_assay

  ## TMB
  if (!is.null(report$content$tmb)) {
    if ("sample_estimate" %in% names(report$content$tmb)) {
      colnames(report$content$tmb$sample_estimate) <-
        toupper(colnames(report$content$tmb$sample_estimate))
      excel_sheets[['TMB']] <- report$content$tmb$sample_estimate
    }
  }

  ## Mutational signatures
  if (!is.null(report$content$mutational_signatures)) {
    if (report$content$mutational_signatures$missing_data == FALSE) {
      if ("result" %in% names(report$content$mutational_signatures)) {
        if ("tsv" %in% names(report$content$mutational_signatures$result)) {
          colnames(report$content$mutational_signatures$result$tsv) <-
            toupper(colnames(report$content$mutational_signatures$result$tsv))

          excel_sheets[['MUTATIONAL_SIGNATURE']] <-
            report$content$mutational_signatures$result$tsv
        }
      }
    }
  }

  ## Kataegis events
  if (!is.null(report$content$kataegis$events)) {
    if (report$content$kataegis$eval == TRUE) {
      if (NROW(report$content$kataegis$events) > 0) {
        colnames(report$content$kataegis$events) <-
          toupper(colnames(report$content$kataegis$events))

        excel_sheets[['KATAEGIS_EVENTS']] <-
          report$content$kataegis$events
      }
    }
  }


  ## MSI
  if (!is.null(report$content$msi)) {
    if (report$content$msi$missing_data == FALSE) {
      if ("prediction" %in% names(report$content$msi)) {
        if ("msi_stats" %in% names(report$content$msi$prediction)) {
          colnames(report$content$msi$prediction$msi_stats) <-
            toupper(colnames(report$content$msi$prediction$msi_stats))

          excel_sheets[['MSI']] <-
            report$content$msi$prediction$msi_stats
        }
      }
    }
  }

  ## Immune contexture
  if (!is.null(report$content$expression) &&
     isTRUE(report$content$expression$eval)) {
    if ("immune_contexture" %in% names(report$content$expression)) {

      colnames(report$content$expression$immune_contexture) <-
        toupper(colnames(report$content$expression$immune_contexture))

      excel_sheets[['RNA_IMMUNE_CONTEXTURE']] <-
        report$content$expression$immune_contexture
    }

    if ("outliers" %in% names(report$content$expression)) {

      excel_sheets[['RNA_EXPRESSION_OUTLIERS']] <-
        report$content$expression$outliers |>
        remove_cols_from_df(
          cnames = c("GENENAME",
                     "CANCERGENE_EVIDENCE",
                     "TARGETED_INHIBITORS_ALL",
                     "ENSEMBL_GENE_ID")) |>
        dplyr::select(-dplyr::any_of(
          c("GENENAME",
            "CANCERGENE_EVIDENCE",
            "TARGETED_INHIBITORS_ALL",
            "ENSEMBL_GENE_ID")))
    }
  }

  ## Copy number alterations
  if (!is.null(report$content$cna) &&
      !is.null(callsets[['cna']]) &&
     isTRUE(report$content$cna$eval)) {

    if (NROW(callsets[['cna']]$variant) > 0) {
      excel_sheets[['SOMATIC_CNA']] <- as.data.frame(
        callsets[['cna']]$variant |>
          dplyr::select(dplyr::any_of(tsv_cols$cna)) |>
          dplyr::select(-dplyr::any_of("BIOMARKER_MATCH")) |>
          dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
          dplyr::arrange(
            .data$ACTIONABILITY_TIER,
            .data$TARGETED_INHIBITORS_ALL2)
      )

      ## Evidence items - biomarkers
      excel_sheets[['SOMATIC_CNA_BIOMARKER']] <- data.frame()

      for (clnsig in names(bm_categories)) {
        if(clnsig == "diagnostic_negative") {
          next
        }
        if(NROW(callsets[['cna']]$bm_evidence[[clnsig]]$eitems) > 0){
          excel_sheets[['SOMATIC_CNA_BIOMARKER']] <-
            dplyr::bind_rows(
              excel_sheets[['SOMATIC_CNA_BIOMARKER']],
              callsets[['cna']]$bm_evidence[[clnsig]]$eitems |>
                dplyr::mutate(
                  TIER = switch(
                    clnsig,
                    therapeutic_sensitivity = paste0("T", .data$ACTIONABILITY_TIER),
                    therapeutic_resistance = paste0("R", .data$ACTIONABILITY_TIER),
                    prognostic_poor = paste0("PP", .data$ACTIONABILITY_TIER),
                    prognostic_better = paste0("PB", .data$ACTIONABILITY_TIER),
                    diagnostic_positive = paste0("D", .data$ACTIONABILITY_TIER),
                    as.character(.data$ACTIONABILITY_TIER)
                  )
                ) |>
                dplyr::select(
                  -dplyr::any_of(c("BM_VARIANT_ID",
                     "ACTIONABILITY_TIER",
                     "BM_REFERENCE",
                     "BM_EVIDENCE_LEVEL_FULL",
                     "BM_EVIDENCE_DIRECTION"))
                ) |>
                dplyr::mutate(BM_MOLECULAR_PROFILE = strip_html(
                  .data$BM_MOLECULAR_PROFILE
                )) |>
                dplyr::mutate(SAMPLE_ID = sample_id) |>
                dplyr::left_join(
                  sample_alteration[['cna']],
                  by = c("VAR_ID",
                         "ENTREZGENE",
                         "VARIANT_CLASS")
                ) |>
                dplyr::select(
                  c("SAMPLE_ID",
                    "VARIANT_CLASS",
                    "VAR_ID",
                    "SAMPLE_ALTERATION",
                    "TIER"),
                  dplyr::everything()
                )
            )
        }
      }


      if (NROW(excel_sheets[['SOMATIC_CNA_BIOMARKER']]) > 0) {
        excel_sheets[['SOMATIC_CNA_BIOMARKER']] <-
          excel_sheets[['SOMATIC_CNA_BIOMARKER']] |>
          dplyr::distinct() |>
          dplyr::mutate(
            .tier_order = dplyr::case_when(
              startsWith(.data$TIER, "T") ~ 1L,
              startsWith(.data$TIER, "R") ~ 2L,
              startsWith(.data$TIER, "D") ~ 3L,
              startsWith(.data$TIER, "P") ~ 4L,
              TRUE ~ 5L),
            BM_ACTIONABILITY_SUPPORT = factor(
              .data$BM_ACTIONABILITY_SUPPORT,
              levels = c("tier-defining", "additional"))) |>
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$.tier_order,
            .data$TIER,
            .data$VAR_ID,
            .data$SAMPLE_ALTERATION,
            .data$BM_ACTIONABILITY_SUPPORT,
            .data$BM_EVIDENCE_LEVEL,
            dplyr::desc(.data$BM_RATING)) |>
          dplyr::mutate(BM_ACTIONABILITY_SUPPORT = as.character(
            .data$BM_ACTIONABILITY_SUPPORT)) |>
          dplyr::select(-.data$.tier_order)
      }
    }
  }


  ## RNA fusions
  if (!is.null(report$content$fusion) &&
      !is.null(callsets[['fusion']]) &&
      isTRUE(report$content$fusion$eval)) {

    if (NROW(callsets[['fusion']]$variant) > 0) {
      excel_sheets[['RNA_FUSION']] <- as.data.frame(
        callsets[['fusion']]$variant |>
          dplyr::select(dplyr::any_of(tsv_cols$fusion)) |>
          dplyr::select(-dplyr::any_of("BIOMARKER_MATCH")) |>
          dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
          dplyr::arrange(
            .data$ACTIONABILITY_TIER)
      )

      ## Evidence items - biomarkers
      excel_sheets[['RNA_FUSION_BIOMARKER']] <- data.frame()

      for (clnsig in names(bm_categories)) {
        if(clnsig == "diagnostic_negative") {
          next
        }
        if(NROW(callsets[['fusion']]$bm_evidence[[clnsig]]$eitems) > 0){
          excel_sheets[['RNA_FUSION_BIOMARKER']] <-
            dplyr::bind_rows(
              excel_sheets[['RNA_FUSION_BIOMARKER']],
              callsets[['fusion']]$bm_evidence[[clnsig]]$eitems |>
                dplyr::mutate(
                  TIER = switch(
                    clnsig,
                    therapeutic_sensitivity = paste0("T", .data$ACTIONABILITY_TIER),
                    therapeutic_resistance = paste0("R", .data$ACTIONABILITY_TIER),
                    prognostic_poor = paste0("PP", .data$ACTIONABILITY_TIER),
                    prognostic_better = paste0("PB", .data$ACTIONABILITY_TIER),
                    diagnostic_positive = paste0("D", .data$ACTIONABILITY_TIER),
                    as.character(.data$ACTIONABILITY_TIER)
                  )
                ) |>
                dplyr::select(
                  -dplyr::any_of(
                    c("BM_VARIANT_ID",
                     "ACTIONABILITY_TIER",
                     "BM_REFERENCE",
                     "BM_EVIDENCE_LEVEL_FULL",
                     "BM_EVIDENCE_DIRECTION"))
                ) |>
                dplyr::mutate(BM_MOLECULAR_PROFILE = strip_html(
                  .data$BM_MOLECULAR_PROFILE
                )) |>
                dplyr::mutate(SAMPLE_ID = sample_id) |>
                dplyr::left_join(
                  ## Fusion biomarker items carry a single-gene ENTREZGENE
                  ## (e.g. "238" for ALK) while the variant table stores the
                  ## resolved pair (e.g. "238::238"). Joining on ENTREZGENE
                  ## fails for single-gene fusion biomarkers and leaves
                  ## SAMPLE_ALTERATION empty. Join on (VAR_ID, VARIANT_CLASS)
                  ## only; drop ENTREZGENE from the right side to avoid
                  ## ENTREZGENE.x / ENTREZGENE.y column collisions.
                  dplyr::select(
                    sample_alteration[['fusion']],
                    -dplyr::any_of("ENTREZGENE")),
                  by = c("VAR_ID",
                         "VARIANT_CLASS")
                ) |>
                dplyr::select(
                  c("SAMPLE_ID",
                    "VARIANT_CLASS",
                    "VAR_ID",
                    "SAMPLE_ALTERATION",
                    "TIER"),
                  dplyr::everything()
                )
            )
        }
      }


      if (NROW(excel_sheets[['RNA_FUSION_BIOMARKER']]) > 0) {
        excel_sheets[['RNA_FUSION_BIOMARKER']] <-
          excel_sheets[['RNA_FUSION_BIOMARKER']] |>
          dplyr::distinct() |>
          dplyr::mutate(
            .tier_order = dplyr::case_when(
              startsWith(.data$TIER, "T") ~ 1L,
              startsWith(.data$TIER, "R") ~ 2L,
              startsWith(.data$TIER, "D") ~ 3L,
              startsWith(.data$TIER, "P") ~ 4L,
              TRUE ~ 5L),
            BM_ACTIONABILITY_SUPPORT = factor(
              .data$BM_ACTIONABILITY_SUPPORT,
              levels = c("tier-defining", "additional"))) |>
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$.tier_order,
            .data$TIER,
            .data$VAR_ID,
            .data$SAMPLE_ALTERATION,
            .data$BM_ACTIONABILITY_SUPPORT,
            .data$BM_EVIDENCE_LEVEL,
            dplyr::desc(.data$BM_RATING)) |>
          dplyr::mutate(BM_ACTIONABILITY_SUPPORT = as.character(
            .data$BM_ACTIONABILITY_SUPPORT)) |>
          dplyr::select(-.data$.tier_order)
      }
    }
  }


  ## SNVs/InDels
  if (!is.null(report$content$snv_indel) &&
      !is.null(callsets$snv_indel) &&
     isTRUE(report$content$snv_indel$eval)) {

    snv_indel_cols <- tsv_cols$snv_indel
    if (report$settings$conf$other$retained_vcf_info_tags != "None") {
      snv_indel_cols <- c(
        snv_indel_cols,
        stringr::str_split(
          report$settings$conf$other$retained_vcf_info_tags, ",")[[1]]
      )
    }


    if (NROW(callsets[['snv_indel']]$variant) > 0) {
      excel_sheets[['SOMATIC_SNV_INDEL']] <-
        callsets[['snv_indel']]$variant |>
        dplyr::select(
          dplyr::any_of(snv_indel_cols)) |>
        dplyr::select(-dplyr::any_of("BIOMARKER_MATCH")) |>
        ## Limit Excel output to exonic variants, as well
        ## as any variant in top tiers of actionability (even if non-exonic)
        dplyr::filter(
          .data$EXONIC_STATUS == "exonic" |
            (.data$EXONIC_STATUS == "nonexonic" &
               !is.na(.data$ACTIONABILITY_TIER) &
               .data$ACTIONABILITY_TIER <= 3)) |>
        dplyr::select(
          -dplyr::any_of(
            c("BIOMARKER_MATCH","VEP_ALL_CSQ")))


      ## Excel workbook won't allow column names that are identical (case insensitive)
      ## - rename HGVSP to HGVSp_short
      if ("HGVSP" %in% colnames(excel_sheets[['SOMATIC_SNV_INDEL']]) &
         "HGVSp" %in% colnames(excel_sheets[['SOMATIC_SNV_INDEL']])) {
        excel_sheets[['SOMATIC_SNV_INDEL']] <-
          excel_sheets[['SOMATIC_SNV_INDEL']] |>
          dplyr::rename(HGVSp_short = "HGVSP")
      }

      ## Evidence items - biomarkers
      excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] <- data.frame()

      for (clnsig in names(bm_categories)) {
        if(clnsig == "diagnostic_negative") {
          next
        }
        if(NROW(callsets[['snv_indel']]$bm_evidence[[clnsig]]$eitems) > 0){
          excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] <-
            dplyr::bind_rows(
              excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']],
              callsets[['snv_indel']]$bm_evidence[[clnsig]]$eitems |>
                dplyr::mutate(
                  TIER = switch(
                    clnsig,
                    therapeutic_sensitivity = paste0("T", .data$ACTIONABILITY_TIER),
                    therapeutic_resistance = paste0("R", .data$ACTIONABILITY_TIER),
                    prognostic_poor = paste0("PP", .data$ACTIONABILITY_TIER),
                    prognostic_better = paste0("PB", .data$ACTIONABILITY_TIER),
                    diagnostic_positive = paste0("D", .data$ACTIONABILITY_TIER),
                    as.character(.data$ACTIONABILITY_TIER)
                  )
                ) |>
                dplyr::select(
                  -dplyr::any_of(c("BM_VARIANT_ID",
                     "ACTIONABILITY_TIER",
                     "BM_REFERENCE",
                     "BM_EVIDENCE_LEVEL_FULL",
                     "BM_EVIDENCE_DIRECTION"))
                ) |>
                dplyr::mutate(BM_MOLECULAR_PROFILE = strip_html(
                  .data$BM_MOLECULAR_PROFILE
                )) |>
                dplyr::mutate(SAMPLE_ID = sample_id) |>
                dplyr::left_join(
                  sample_alteration[['snv_indel']],
                  by = c("VAR_ID","ENTREZGENE","VARIANT_CLASS")
                ) |>
                dplyr::select(
                  c("SAMPLE_ID",
                    "VARIANT_CLASS",
                    "VAR_ID",
                    "SAMPLE_ALTERATION",
                    "TIER"),
                  dplyr::everything()
                )
            )
        }
      }


      if (NROW(excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']]) > 0) {
        excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] <-
          excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] |>
          dplyr::distinct() |>
          dplyr::mutate(
            .tier_order = dplyr::case_when(
              startsWith(.data$TIER, "T") ~ 1L,
              startsWith(.data$TIER, "R") ~ 2L,
              startsWith(.data$TIER, "D") ~ 3L,
              startsWith(.data$TIER, "P") ~ 4L,
              TRUE ~ 5L),
            BM_ACTIONABILITY_SUPPORT = factor(
              .data$BM_ACTIONABILITY_SUPPORT,
              levels = c("tier-defining", "additional"))) |>
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$.tier_order,
            .data$TIER,
            .data$VAR_ID,
            .data$SAMPLE_ALTERATION,
            .data$BM_ACTIONABILITY_SUPPORT,
            .data$BM_EVIDENCE_LEVEL,
            dplyr::desc(.data$BM_RATING)) |>
          dplyr::mutate(BM_ACTIONABILITY_SUPPORT = as.character(
            .data$BM_ACTIONABILITY_SUPPORT)) |>
          dplyr::select(-.data$.tier_order)
      }
    }
  }

  for (e in names(excel_sheets)) {
    if (NROW(excel_sheets[[e]]) == 0) {
      excel_sheets[[e]] <- NULL
    }

    if (NROW(excel_sheets[[e]]) > 30000) {
      log4r_warn(paste0(
        "Excel sheet '", e, "' has more than 30,000 rows - ",
        "only first 30,000 rows will be written to Excel workbook - consider TSV for full data"))
      excel_sheets[[e]] <- utils::head(
        excel_sheets[[e]], 30000)
    }

  }


  return(excel_sheets)

}
