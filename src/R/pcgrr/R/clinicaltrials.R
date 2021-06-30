#' Function that retrieves relevant (interventional based on molecular target)
#' clinical trials for a given tumor type
#' @param pcgr_data PCGR data bundle object
#' @param config PCGR run configurations
#' @param sample_name sample name
#'
#' @return pcg_report_trials data frame with all report elements
#'

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
    pcgr_data[["clinicaltrials"]][["trials"]] %>%
    dplyr::filter(primary_site ==
                    config[["t_props"]][["tumor_type"]])

  if (nrow(pcg_report_trials[["trials"]]) > 0) {

    pcg_report_trials[["trials"]] <- pcg_report_trials[["trials"]] %>%
      dplyr::select(nct_id, title, overall_status,
                    cui_link, intervention_link,
                    phase, start_date,
                    primary_completion_date, cui_name,
                    intervention,
                    intervention_target,
                    biomarker_context,
                    chromosome_abnormality,
                    clinical_context,
                    world_region,
                    metastases, gender,
                    minimum_age, maximum_age, phase,
                    n_primary_cancer_sites,
                    study_design_primary_purpose) %>%
      dplyr::rename(condition_raw = cui_name,
                    condition = cui_link,
                    intervention2 = intervention_link,
                    intervention_raw = intervention,
                    biomarker_index = biomarker_context,
                    keyword = clinical_context,
                    chrom_abnormalities = chromosome_abnormality,
                    metastases_index = metastases) %>%
      dplyr::rename(intervention = intervention2) %>%
      magrittr::set_colnames(toupper(names(.))) %>%
      dplyr::arrange(N_PRIMARY_CANCER_SITES,
                     OVERALL_STATUS, desc(START_DATE),
                     desc(nchar(BIOMARKER_INDEX),
                          STUDY_DESIGN_PRIMARY_PURPOSE)) %>%
      dplyr::select(-c(N_PRIMARY_CANCER_SITES, STUDY_DESIGN_PRIMARY_PURPOSE))

    if (nrow(pcg_report_trials[["trials"]]) > 2000) {
      pcg_report_trials[["trials"]] <-
        head(pcg_report_trials[["trials"]], 2000)
    }

  }else{
    pcg_report_trials[["missing_data"]] <- T
  }
  return(pcg_report_trials)
}
