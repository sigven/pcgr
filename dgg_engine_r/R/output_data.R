#' Function that produces the contents of sheets for an Excel report
#' of PCGR output
#'
#' @param report PCGR report object
#'
#' @export
#'
get_excel_sheets <- function(report = NULL){

  if(is.null(report)){
    stop("report must be provided")
  }

  if(!((base::typeof(report$content) == "list") == TRUE)){
    stop("report$content must be a list")
  }

  ## Sample and assay properties
  excel_sheets <- list()
  excel_sheets[['SAMPLE_ASSAY']] <- data.frame()
  sample_id <- report$settings$sample_id

  if(!((base::typeof(report$content$sample_properties) == "list") == TRUE)){
    stop("report$content$sample_properties must be a list")
  }
  init_sample_assay <-
    dplyr::bind_rows(
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'SAMPLE',
        PROPERTY = 'SITE', VALUE = NA
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

  for(elem in c('SITE','TUMOR_PURITY','TUMOR_PLOIDY')){
    if(tolower(elem) %in% names(report$content$sample_properties)){
        init_sample_assay[init_sample_assay$PROPERTY == elem, 'VALUE'] <-
          report$content$sample_properties[[tolower(elem)]]
      }
    }

  for(elem in c('TYPE','MODE','EFFECTIVE_TARGET_SIZE_MB')){
    if(tolower(elem) %in% names(report$content$assay_properties)){
      init_sample_assay[init_sample_assay$PROPERTY == elem, 'VALUE'] <-
        report$content$assay_properties[[tolower(elem)]]
    }
  }

  excel_sheets[['SAMPLE_ASSAY']] <- init_sample_assay

  ## TMB
  if(!is.null(report$content$tmb)){
    if("sample_estimate" %in% names(report$content$tmb)){
      colnames(report$content$tmb$sample_estimate) <-
        toupper(colnames(report$content$tmb$sample_estimate))
      excel_sheets[['TMB']] <- report$content$tmb$sample_estimate
    }
  }

  ## Mutational signatures
  if(!is.null(report$content$mutational_signatures)){
    if(report$content$mutational_signatures$missing_data == FALSE){
      if("result" %in% names(report$content$mutational_signatures)){
        if("tsv" %in% names(report$content$mutational_signatures$result)){
          colnames(report$content$mutational_signatures$result$tsv) <-
            toupper(colnames(report$content$mutational_signatures$result$tsv))

          excel_sheets[['MUTATIONAL_SIGNATURE']] <-
            report$content$mutational_signatures$result$tsv
        }
      }
    }
  }

  ## Kataegis events
  if(!is.null(report$content$kataegis$events)){
    if(report$content$kataegis$eval == TRUE){
      if(NROW(report$content$kataegis$events) > 0){
        colnames(report$content$kataegis$events) <-
          toupper(colnames(report$content$kataegis$events))

        excel_sheets[['KATAEGIS_EVENTS']] <-
          report$content$kataegis$events
      }
    }
  }


  ## MSI
  if(!is.null(report$content$msi)){
    if(report$content$msi$missing_data == FALSE){
      if("prediction" %in% names(report$content$msi)){
        if("msi_stats" %in% names(report$content$msi$prediction)){
          colnames(report$content$msi$prediction$msi_stats) <-
            toupper(colnames(report$content$msi$prediction$msi_stats))

          excel_sheets[['MSI']] <-
            report$content$msi$prediction$msi_stats
        }
      }
    }
  }

  ## Immune contexture
  if(!is.null(report$content$expression) &
     report$content$expression$eval == TRUE){
    if("immune_contexture" %in% names(report$content$expression)){

      colnames(report$content$expression$immune_contexture) <-
        toupper(colnames(report$content$expression$immune_contexture))

      excel_sheets[['RNA_IMMUNE_CONTEXTURE']] <-
        report$content$expression$immune_contexture
    }

    if("outliers" %in% names(report$content$expression)){

      excel_sheets[['RNA_EXPRESSION_OUTLIERS']] <-
        report$content$expression$outliers |>
        dplyr::select(-dplyr::any_of(
          c("GENENAME","CANCERGENE_EVIDENCE",
            "TARGETED_INHIBITORS_ALL","ENSEMBL_GENE_ID")))
    }
  }

  ## Copy number alterations
  if(!is.null(report$content$cna) &
     report$content$cna$eval == TRUE){

    if(NROW(report$content$cna$callset$variant) > 0){
      excel_sheets[['SOMATIC_CNA']] <- as.data.frame(
        report$content$cna$callset$variant |>
          dplyr::select(dplyr::any_of(pcgrr::tsv_cols$cna)) |>
          dplyr::select(-dplyr::any_of("BIOMARKER_MATCH")) |>
          dplyr::filter(!is.na(.data$ACTIONABILITY_TIER))
      )

      ## Evidence items - biomarkers
      excel_sheets[['SOMATIC_CNA_BIOMARKER']] <- data.frame()
      i <- 1
      while(i <= 2){
        tier_data <-
          get_dt_tables(
            rep = report,
            tier = i,
            variant_class = "cna")
        if(NROW(tier_data$by_eitem) > 0){
          edata <- tier_data$by_eitem |>
            dplyr::mutate(BM_REFERENCE = strip_html(
              .data$BM_REFERENCE
            )) |>
            dplyr::mutate(BM_MOLECULAR_PROFILE = strip_html(
              .data$BM_MOLECULAR_PROFILE
            )) |>
            dplyr::select(-c("BM_CONTEXT")) |>
            dplyr::mutate(SAMPLE_ID = sample_id,
                          ACTIONABILITY_TIER = i) |>
            dplyr::rename(SAMPLE_ALTERATION = "MOLECULAR_ALTERATION") |>
            dplyr::group_by(dplyr::across(-c("BM_PRIMARY_SITE"))) |>
            dplyr::summarise(
              BM_PRIMARY_SITE = paste(
                unique(.data$BM_PRIMARY_SITE), collapse = ", "),
              .groups = "drop") |>
            dplyr::select(
              c("SAMPLE_ID","SAMPLE_ALTERATION", "ACTIONABILITY_TIER"),
              dplyr::everything()
            )
          excel_sheets[['SOMATIC_CNA_BIOMARKER']] <- dplyr::bind_rows(
            excel_sheets[['SOMATIC_CNA_BIOMARKER']], edata)
        }
        i <- i + 1
      }
      if(NROW(excel_sheets[['SOMATIC_CNA_BIOMARKER']]) > 0){
        excel_sheets[['SOMATIC_CNA_BIOMARKER']] <-
          excel_sheets[['SOMATIC_CNA_BIOMARKER']] |>
          dplyr::distinct() |>
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$ACTIONABILITY_TIER,
            .data$BM_EVIDENCE_LEVEL,
            dplyr::desc(.data$BM_RATING))
      }
    }
  }

  ## SNVs/InDels
  if(!is.null(report$content$snv_indel) &
     report$content$snv_indel$eval == TRUE){

    snv_indel_cols <- pcgrr::tsv_cols$snv_indel
    if(report$settings$conf$other$retained_vcf_info_tags != "None"){
      snv_indel_cols <- c(
        snv_indel_cols,
        stringr::str_split(
          report$settings$conf$other$retained_vcf_info_tags, ",")[[1]]
      )
    }

    if(NROW(report$content$snv_indel$callset$variant) > 0){
      excel_sheets[['SOMATIC_SNV_INDEL']] <-
        report$content$snv_indel$callset$variant |>
        dplyr::select(
          dplyr::any_of(snv_indel_cols)) |>
        ## Limit Excel output to exonic variants, as well
        ## as any actionable variant (TIER I-II)
        dplyr::filter(.data$EXONIC_STATUS == "exonic" |
                        (.data$EXONIC_STATUS == "nonexonic" &
                           !is.na(.data$ACTIONABILITY_TIER) &
                           .data$ACTIONABILITY_TIER <= 2)) |>
        dplyr::select(
          -dplyr::any_of(c("BIOMARKER_MATCH","VEP_ALL_CSQ")))

      ## Evidence items - biomarkers
      excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] <- data.frame()
      i <- 1
      while(i <= 2){
        tier_data <-
          get_dt_tables(
            rep = report,
            tier = i,
            variant_class = "snv_indel")
        if(NROW(tier_data$by_eitem) > 0){
          edata <- tier_data$by_eitem |>
            dplyr::mutate(BM_REFERENCE = strip_html(
              .data$BM_REFERENCE
            )) |>
            dplyr::mutate(BM_MOLECULAR_PROFILE = strip_html(
              .data$BM_MOLECULAR_PROFILE
            )) |>
            dplyr::select(-c("BM_CONTEXT")) |>
            dplyr::mutate(SAMPLE_ID = sample_id,
                          ACTIONABILITY_TIER = i) |>
            dplyr::rename(SAMPLE_ALTERATION = "MOLECULAR_ALTERATION") |>
            dplyr::group_by(dplyr::across(-c("BM_PRIMARY_SITE"))) |>
            dplyr::summarise(
              BM_PRIMARY_SITE = paste(
                unique(.data$BM_PRIMARY_SITE), collapse = ", "),
              .groups = "drop") |>
            dplyr::select(
              c("SAMPLE_ID","SAMPLE_ALTERATION", "ACTIONABILITY_TIER"),
              dplyr::everything()
            )

          excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] <- dplyr::bind_rows(
            excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']], edata)
        }
        i <- i + 1
      }
      if(NROW(excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']]) > 0){
        excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] <-
          excel_sheets[['SOMATIC_SNV_INDEL_BIOMARKER']] |>
          dplyr::distinct() |>
          dplyr::arrange(
            .data$SAMPLE_ID,
            .data$ACTIONABILITY_TIER,
            .data$BM_EVIDENCE_LEVEL,
            dplyr::desc(.data$BM_RATING))
      }
    }
  }

  for(e in names(excel_sheets)){
    if(NROW(excel_sheets[[e]]) == 0){
      excel_sheets[[e]] <- NULL
    }
  }


  return(excel_sheets)

}
