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

  ## Sample and assay properties
  excel_sheets <- list()
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
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$TIER,
            .data$BM_EVIDENCE_LEVEL,
            .data$BM_PRIMARY_SITE,
            dplyr::desc(.data$BM_RATING))
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
                  sample_alteration[['fusion']],
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


      if (NROW(excel_sheets[['RNA_FUSION_BIOMARKER']]) > 0) {
        excel_sheets[['RNA_FUSION_BIOMARKER']] <-
          excel_sheets[['RNA_FUSION_BIOMARKER']] |>
          dplyr::distinct() |>
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$TIER,
            .data$BM_EVIDENCE_LEVEL,
            .data$BM_PRIMARY_SITE,
            dplyr::desc(.data$BM_RATING))
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
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$TIER,
            .data$BM_EVIDENCE_LEVEL,
            .data$BM_PRIMARY_SITE,
            dplyr::desc(.data$BM_RATING))
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
