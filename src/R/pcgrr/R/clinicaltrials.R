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


  pcg_report_trials <- pcgrr::init_report(config, sample_name,
                                          class = "clinicaltrials")
  pcg_report_trials[["eval"]] <- T

  pcg_report_trials[["trials"]] <-
    pcgr_data[["clinicaltrials"]][["trials"]] %>%
    dplyr::filter(primary_site ==
                    config[["t_props"]][["tumor_type"]])

  if (nrow(pcg_report_trials[["trials"]]) > 0) {
    pcg_report_trials[["trials"]] <- pcg_report_trials[["trials"]] %>%
      dplyr::filter(
        (!is.na(keyword) |
           intervention_n_target > 0 |
           !is.na(biomarker)) &
          (stringr::str_detect(
            study_design_primary_purpose,
            "Treatment|Other|Prevention|Diagnostic|Basic Science")) &
          (stringr::str_detect(
            overall_status,
            "ecruiting|Enrolling|Unknown|Suspended|Withdrawn|Completed")))

    if (nrow(pcg_report_trials[["trials"]]) > 0) {

      pcg_report_trials[["trials"]] <- pcg_report_trials[["trials"]] %>%
        dplyr::select(nct_id, title, overall_status,
                      cui_link, keyword, intervention_link,
                      phase, start_date,
                      primary_completion_date, cui_name,
                      intervention, intervention_target,
                      biomarker, biomarker_support,
                      metastases, gender,
                      minimum_age, maximum_age, phase,
                      num_primary_sites,
                      study_design_primary_purpose) %>%
        dplyr::rename(condition_raw = cui_name,
                      condition = cui_link,
                      intervention2 = intervention_link,
                      intervention_raw = intervention,
                      biomarker_index = biomarker,
                      biomarker_index_support = biomarker_support,
                      metastases_index = metastases) %>%
        dplyr::rename(intervention = intervention2) %>%
        magrittr::set_colnames(toupper(names(.))) %>%
        dplyr::select(-BIOMARKER_INDEX_SUPPORT) %>%
        dplyr::arrange(NUM_PRIMARY_SITES,
                       OVERALL_STATUS, desc(START_DATE),
                       desc(nchar(BIOMARKER_INDEX),
                            STUDY_DESIGN_PRIMARY_PURPOSE)) %>%
        dplyr::select(-c(NUM_PRIMARY_SITES, STUDY_DESIGN_PRIMARY_PURPOSE))

      if (nrow(pcg_report_trials[["trials"]]) > 2000) {
        pcg_report_trials[["trials"]] <-
          head(pcg_report_trials[["trials"]], 2000)
      }


    }else{
      pcg_report_trials[["missing_data"]] <- T
    }

  }else{
    pcg_report_trials[["missing_data"]] <- T
  }
  return(pcg_report_trials)
}
