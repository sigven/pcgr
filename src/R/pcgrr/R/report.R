
#' Function that initates PCGR report object
#'
#' @param pcgr_config Object with PCGR configuration parameters
#' @param cna_segments_tsv name of CNA segments file (tab-separated values)
#' @param sample_name sample identifier
#' @param pcgr_version PCGR software version
#' @param genome_assembly human genome assembly version
#' @param class report analysis section (NULL defaults to full report)
#' @param pcgr_data pcgr data object

init_pcg_report <- function(pcgr_config = NULL, sample_name = 'SampleX', pcgr_version = '0.6.0', genome_assembly = 'grch37', class = NULL, pcgr_data = NULL, type = 'somatic'){


  pcg_report <- list()
  pcg_report[['pcgr_config']] <- pcgr_config
  pcg_report[['sample_name']] <- sample_name
  pcg_report[['genome_assembly']] <- genome_assembly
  pcg_report[['pcgr_version']] <- pcgr_version
  pcg_report[['pcgr_db_release']] <- pcgr_data$pcgr_db_release

  if(type == 'predispose'){
    analysis_element <- 'snv_indel'
    pcg_report[[analysis_element]] <- list()
    pcg_report[[analysis_element]][['eval']] <- FALSE
    pcg_report[[analysis_element]][['variant_display']] <- list()
    pcg_report[[analysis_element]][['variant_set']] <- list()
    pcg_report[[analysis_element]][['variant_statistic']] <- list()
    pcg_report[[analysis_element]][['zero']] <- FALSE

    cancer_genes <- pcgrr::list_to_df(pcgr_config$cancer_predisposition_genes) %>% dplyr::filter(list.element == T) %>% dplyr::select(name) %>% dplyr::rename(symbol = name)
    pcg_report[[analysis_element]][['predisposition_genes']] <- cancer_genes

    for(t in c('tier1','tier2','tier3A','tier3B','gwas')){
      pcg_report[[analysis_element]][['variant_display']][[t]] <- data.frame()
      if(t != 'tier3B' & t != 'gwas'){
        pcg_report[[analysis_element]][['variant_display']][[t]] <- list()
        for(c in c('cancer_phenotype','noncancer_phenotype')){
          pcg_report[[analysis_element]][['variant_display']][[t]][[c]] <- data.frame()
        }
      }
      pcg_report[[analysis_element]][['variant_set']][[t]] <- data.frame()
    }
    pcg_report[[analysis_element]][['variant_set']][['tsv']] <- data.frame()

    for(t in c('n','n_snv','n_indel','n_coding','n_noncoding')){
      pcg_report[[analysis_element]][['variant_statistic']][[t]] <- 0
    }

    if(!is.null(pcg_report[['pcgr_config']][['popgen']])){
      if(pcg_report[['pcgr_config']][['popgen']][['pop_tgp']] != ""){
        pop_tag_info <- pcgrr::get_population_tag(pcg_report[['pcgr_config']][['popgen']][['pop_tgp']], db = "1KG")
        pcg_report[['pcgr_config']][['popgen']][['vcftag_tgp']] <- pop_tag_info$vcf_tag
        pcg_report[['pcgr_config']][['popgen']][['popdesc_tgp']] <- pop_tag_info$pop_description
      }
      if(pcgr_config[['popgen']][['pop_gnomad']] != ""){
        pop_tag_info <- pcgrr::get_population_tag(pcgr_config[['popgen']][['pop_gnomad']], db = "GNOMAD")
        pcg_report[['pcgr_config']][['popgen']][['vcftag_gnomad']] <- pop_tag_info$vcf_tag
        pcg_report[['pcgr_config']][['popgen']][['popdesc_gnomad']] <- pop_tag_info$pop_description
      }
    }

    pcg_report[['summary']] <- list()
    for(t in c("tier1","tier2","tier3A","tier3B","gwas")){
      pcg_report[['summary']][[t]] <- data.frame()
    }
    return(pcg_report)
  }

  pcg_report[['tier_model']] <- pcgr_config$tier_model$tier_model

  pcg_report[['tumor_class']] <- 'Not defined'
  tumor_types_set <- pcgrr::list_to_df(pcgr_config$tumor_type) %>% dplyr::filter(list.element == T) %>% dplyr::select(name)
  if(nrow(tumor_types_set) > 0){
    pcg_report[['tumor_class']] <- paste0(tumor_types_set$name,collapse=", ")
  }
  for(analysis_element in c('snv_indel','tmb','msi','cna','cna_plot','m_signature','sequencing_mode','tumor_only','value_box')){
    pcg_report[[analysis_element]] <- list()
    pcg_report[[analysis_element]][['eval']] <- FALSE

    if(analysis_element == 'sequencing_mode'){
      pcg_report[[analysis_element]][['tumor_only']] <- FALSE
      pcg_report[[analysis_element]][['mode']] <- 'Tumor vs. control'
    }

    if(analysis_element == 'value_box'){
      pcg_report[[analysis_element]][['tmb_tertile']] <- 'TMB:\nNot determined'
      pcg_report[[analysis_element]][['msi']] <- 'MSI:\nNot determined'
      pcg_report[[analysis_element]][['scna']] <- 'SCNA:\nNot determined'
      pcg_report[[analysis_element]][['tier1']] <- 'Tier 1 variants:\nNot determined'
      pcg_report[[analysis_element]][['tier2']] <- 'Tier 2 variants:\nNot determined'
      pcg_report[[analysis_element]][['signatures']] <- 'Mutational signatures:\nNot determined'
    }

    if(analysis_element == 'snv_indel' | analysis_element == 'cna'){
      pcg_report[[analysis_element]][['clinical_evidence_item']] <- list()
      pcg_report[[analysis_element]][['variant_display']] <- list()
      pcg_report[[analysis_element]][['variant_set']] <- list()
      pcg_report[[analysis_element]][['variant_statistic']] <- list()
      pcg_report[[analysis_element]][['zero']] <- FALSE
      for(tumorclass in c('any_tumortype','other_tumortype','specific_tumortype')){
        pcg_report[[analysis_element]][['clinical_evidence_item']][[tumorclass]] <- list()
        for(evidence_type in c('prognostic','diagnostic','predictive')){
          for(evidence_level in c('A_B','C_D_E','any')){
            pcg_report[[analysis_element]][['clinical_evidence_item']][[tumorclass]][[evidence_type]][[evidence_level]] <- data.frame()
          }
        }
      }
      if(analysis_element == 'snv_indel'){
        for(t in c('tier1','tier2','tier3','tier4','noncoding')){
          pcg_report[[analysis_element]][['variant_display']][[t]] <- data.frame()
          if(t == 'tier2' & pcgr_config$tier_model$tier_model == 'pcgr'){
            pcg_report[[analysis_element]][['variant_display']][[t]] <- list()
            for(c in c('hotspot','curated_mutation','predicted_driver')){
              pcg_report[[analysis_element]][['variant_display']][[t]][[c]] <- data.frame()
            }
          }
          if(t == 'tier3' & pcgr_config$tier_model$tier_model == 'pcgr_acmg'){
            pcg_report[[analysis_element]][['variant_display']][[t]] <- list()
            for(c in c('proto_oncogene','tumor_suppressor')){
              pcg_report[[analysis_element]][['variant_display']][[t]][[c]] <- data.frame()
            }
          }
        }
        for(t in c('tier1','tier2','tier3','tier4','noncoding','tsv','tsv_unfiltered','maf','coding','all')){
          pcg_report[[analysis_element]][['variant_set']][[t]] <- data.frame()
        }
        for(t in c('n','n_snv','n_indel','n_coding','n_noncoding','n_tier1','n_tier2','n_tier3','n_tier4')){
          pcg_report[[analysis_element]][['variant_statistic']][[t]] <- 0
        }
      }

      if(analysis_element == 'cna'){
        pcg_report[[analysis_element]][['variant_set']][['cna_print']] <- data.frame()
        for(t in c('n_cna_loss','n_cna_gain')){
          pcg_report[[analysis_element]][['variant_statistic']][[t]] <- data.frame()
        }
        for(t in c('segment','oncogene_gain','tsgene_loss','biomarker')){
          pcg_report[[analysis_element]][['variant_display']][[t]] <- data.frame()
        }
        if(pcgr_config$tier_model$tier_model == 'pcgr_acmg'){
          pcg_report[[analysis_element]][['variant_display']][['biomarkers_tier1']] <- data.frame()
          pcg_report[[analysis_element]][['variant_display']][['biomarkers_tier2']] <- data.frame()
        }
      }
    }
    if(analysis_element == 'm_signature'){
      pcg_report[[analysis_element]][['variant_set']] <- list()
      pcg_report[[analysis_element]][['variant_set']][['all']] <- data.frame()
      pcg_report[[analysis_element]][['missing_data']] <- FALSE
      pcg_report[[analysis_element]][['result']] <- list()
    }
    if(analysis_element == 'tmb'){
      pcg_report[[analysis_element]][['variant_statistic']] <- list()
      pcg_report[[analysis_element]][['variant_statistic']][['n_tmb']] <- 0
      pcg_report[[analysis_element]][['variant_statistic']][['tmb_estimate']] <- 0
      pcg_report[[analysis_element]][['variant_statistic']][['target_size_mb']] <- pcgr_config[['tmb']][['target_size_mb']]
      pcg_report[[analysis_element]][['variant_statistic']][['tmb_tertile']] <- 'TMB - not determined'
    }
    if(analysis_element == 'msi'){
      pcg_report[[analysis_element]][['missing_data']] <- FALSE
      pcg_report[[analysis_element]][['prediction']] <- list()
    }
    if(analysis_element == 'tumor_only'){
      pcg_report[[analysis_element]][['variant_set']] <- list()
      pcg_report[[analysis_element]][['variant_set']][['unfiltered']] <- data.frame()
      pcg_report[[analysis_element]][['variant_set']][['filtered']] <- data.frame()
      pcg_report[[analysis_element]][['variant_statistic']] <- list()

      for(successive_filter in c('unfiltered_n','onekg_n_remain','gnomad_n_remain','dbsnp_n_remain',
                                 'noncoding_n_remain','onekg_frac_remain','gnomad_frac_remain',
                                 'dbsnp_frac_remain','noncoding_frac_remain')){
        pcg_report[[analysis_element]][['variant_statistic']][[successive_filter]] <- 0
      }
    }
  }
  if(!is.null(class)){
    if(!is.null(pcg_report[[class]])){
      return(pcg_report[[class]])
    }
  }

  return(pcg_report)
}

#' Function that initates PCGR report object
#'
#' @param pcg_report PCGR final report
#' @param report_data Object with PCGR report data
#' @param analysis_element section of PCGR report

update_pcg_report <- function(pcg_report, report_data, analysis_element = 'snv_indel'){

  if(!is.null(report_data) & !is.null(pcg_report[[analysis_element]])){
    for(report_elem in names(report_data)){
      if(!is.null(report_data[[report_elem]]) & !is.null(pcg_report[[analysis_element]][[report_elem]])){
         pcg_report[[analysis_element]][[report_elem]] <- report_data[[report_elem]]
      }
    }
  }
  return(pcg_report)
}
