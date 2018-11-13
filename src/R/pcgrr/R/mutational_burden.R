#' Function that tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
#' @return pcg_report_data data frame with all report elements
#'
generate_report_data_tmb <- function(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly = 'hg19'){

  tmb_consequence_pattern <- "^(stop_|start_lost|frameshift_|missense_|synonymous_|inframe_)"
  rlogging::message('------')
  rlogging::message(paste0("Calculating tumor mutational burden"))

  pcg_report_tmb <- pcgrr::init_pcg_report(pcgr_config, sample_name = sample_name, pcgr_version = pcgr_version, genome_assembly = genome_assembly, class = 'tmb')

  pcg_report_tmb[['eval']] <- TRUE
  pcg_report_tmb[['variant_statistic']][['n_tmb']] <- sample_calls %>% dplyr::filter(stringr::str_detect(CONSEQUENCE,tmb_consequence_pattern)) %>% nrow()
  if(pcg_report_tmb[['variant_statistic']][['n_tmb']] > 0 & pcgr_config$mutational_burden$target_size_mb > 0){
    pcg_report_tmb[['variant_statistic']][['tmb_estimate']] <- round(as.numeric(pcg_report_tmb[['variant_statistic']][['n_tmb']] / pcgr_config$mutational_burden$target_size_mb),digits = 2)
    pcg_report_tmb[['variant_statistic']][['target_size_mb']] <- pcgr_config$mutational_burden$target_size_mb
    if(pcg_report_tmb[['variant_statistic']][['tmb_estimate']] <= pcgr_config$mutational_burden$tmb_low_limit){
      pcg_report_tmb[['variant_statistic']][['tmb_tertile']] <- paste0('TMB - Low\n(0 - ',pcgr_config$mutational_burden$tmb_low_limit,' mutations/Mb)')
    }
    if(pcg_report_tmb[['variant_statistic']][['tmb_estimate']] > pcgr_config$mutational_burden$tmb_low_limit & pcg_report_tmb[['variant_statistic']][['tmb_estimate']] <= pcgr_config$mutational_burden$tmb_intermediate_limit){
      pcg_report_tmb[['variant_statistic']][['tmb_tertile']] <- paste0('TMB - Intermediate\n(',pcgr_config$mutational_burden$tmb_low_limit,' - ',pcgr_config$mutational_burden$tmb_intermediate_limit,' mutations/Mb)')
    }
    if(pcg_report_tmb[['variant_statistic']][['tmb_estimate']] > pcgr_config$mutational_burden$tmb_intermediate_limit){
      pcg_report_tmb[['variant_statistic']][['tmb_tertile']] <- paste0('TMB - High\n( > ',pcgr_config$mutational_burden$tmb_intermediate_limit,' mutations/Mb)')
    }
  }
  rlogging::message(paste0("Number of variants for mutational burden analysis: ",pcg_report_tmb[['variant_statistic']][['n_tmb']]))
  rlogging::message(paste0("Estimated mutational burden: ", pcg_report_tmb[['variant_statistic']][['tmb_estimate']]," mutations/Mb"))
  rlogging::message(paste0("Mutational burden tertile: ",stringr::str_replace(pcg_report_tmb[['variant_statistic']][['tmb_tertile']],"\n"," ")))
  rlogging::message('------')

  return(pcg_report_tmb)
}

