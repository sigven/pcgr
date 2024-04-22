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
        PROPERTY = 'PURITY', VALUE = NA
      ),
      data.frame(
        SAMPLE_ID = sample_id, CATEGORY = 'SAMPLE',
        PROPERTY = 'PLOIDY', VALUE = NA
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

  for(elem in c('SITE','PURITY','PLOIDY')){
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

          excel_sheets[['MUTATIONAL_SIGNATURES']] <-
            report$content$mutational_signatures$result$tsv
        }
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

      excel_sheets[['IMMUNE_CONTEXTURE']] <-
        report$content$expression$immune_contexture
    }
  }

  ## Copy number alterations
  if(!is.null(report$content$cna) &
     report$content$cna$eval == TRUE){

    excel_sheets[['CNA']] <- as.data.frame(
      report$content$cna$callset$variant |>
        dplyr::select(dplyr::any_of(pcgrr::tsv_cols$cna)) |>
        dplyr::filter(!is.na(.data$ACTIONABILITY_TIER))
    )
  }

  ## SNVs/InDels
  if(!is.null(report$content$snv_indel) &
     report$content$snv_indel$eval == TRUE){

    snv_indel_cols <- pcgrr::tsv_cols$snv_indel
    if(report$settings$conf$other$retained_vcf_info_tags != "None"){
      snv_indel_cols <- c(
        snv_indel_cols, report$settings$conf$other$retained_vcf_info_tags)
    }

    excel_sheets[['SNV_INDEL']] <- report$content$snv_indel$callset$variant |>
      dplyr::select(
        dplyr::any_of(snv_indel_cols)) |>
      dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
      dplyr::select(-c(dplyr::any_of("VEP_ALL_CSQ")))
  }


  return(excel_sheets)

}
