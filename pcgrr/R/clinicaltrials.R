#' Function that retrieves relevant (interventional based on molecular target)
#' clinical trials for a given tumor type
#' @param pcgr_data PCGR data bundle object
#' @param config PCGR run configurations
#' @param sample_name sample name
#'
#' @return pcg_report_trials data frame with all report elements
#' @export

generate_report_data_trials <- function(pcgr_data, config, sample_name) {


  invisible(assertthat::assert_that(!is.null(pcgr_data)))
  invisible(assertthat::assert_that(!is.null(pcgr_data$clinicaltrials$trials)))
  invisible(assertthat::assert_that(is.data.frame(
    pcgr_data$clinicaltrials$trials)))
  invisible(assertthat::assert_that(
    NROW(pcgr_data$clinicaltrials$trials) > 0))


  pcg_report_trials <- pcgrr::init_report(config = config,
                                          class = "clinicaltrials")
  pcg_report_trials[["eval"]] <- T

  pcg_report_trials[["trials"]] <-
    pcgr_data[["clinicaltrials"]][["trials"]] |>
    dplyr::filter(.data$primary_site ==
                    config[["t_props"]][["tumor_type"]])

  if (nrow(pcg_report_trials[["trials"]]) > 0) {

    pcg_report_trials[["trials"]] <- pcg_report_trials[["trials"]] |>
      dplyr::select(.data$nct_id, .data$title, .data$overall_status,
                    .data$cui_link, .data$intervention_link,
                    .data$phase, .data$start_date,
                    .data$primary_completion_date, .data$cui_name,
                    .data$intervention,
                    .data$intervention_target,
                    .data$biomarker_context,
                    .data$chromosome_abnormality,
                    .data$clinical_context,
                    .data$world_region,
                    .data$metastases, .data$gender,
                    .data$minimum_age, .data$maximum_age, .data$phase,
                    .data$n_primary_cancer_sites,
                    .data$study_design_primary_purpose) |>
      dplyr::rename(condition_raw = .data$cui_name,
                    condition = .data$cui_link,
                    intervention2 = .data$intervention_link,
                    intervention_raw = .data$intervention,
                    biomarker_index = .data$biomarker_context,
                    keyword = .data$clinical_context,
                    chrom_abnormalities = .data$chromosome_abnormality,
                    metastases_index = .data$metastases) |>
      dplyr::rename(intervention = .data$intervention2) |>
      magrittr::set_colnames(toupper(names(.))) |>
      dplyr::arrange(.data$N_PRIMARY_CANCER_SITES,
                     .data$OVERALL_STATUS,
                     dplyr::desc(.data$START_DATE),
                     dplyr::desc(nchar(.data$BIOMARKER_INDEX)),
                     dplyr::desc(.data$STUDY_DESIGN_PRIMARY_PURPOSE)) |>
      dplyr::select(-c(.data$N_PRIMARY_CANCER_SITES, .data$STUDY_DESIGN_PRIMARY_PURPOSE))

    if (nrow(pcg_report_trials[["trials"]]) > 2000) {
      pcg_report_trials[["trials"]] <-
        utils::head(pcg_report_trials[["trials"]], 2000)
    }

  }else{
    pcg_report_trials[["missing_data"]] <- T
  }
  return(pcg_report_trials)
}
