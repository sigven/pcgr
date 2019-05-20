
#' Function that initates PCGR report object
#'
#' @param config Object with configuration parameters
#' @param sample_name sample identifier
#' @param class report analysis section (NULL defaults to full report)
#' @param pcgr_data PCGR software version
#' @param type somatic or predisposition
#' @param virtual_panel_id identifier for virtual panel id

init_pcg_report <- function(config = NULL, sample_name = 'SampleX', class = NULL, pcgr_data = NULL, type = 'somatic', virtual_panel_id = -1){

  report <- list()
  for(elem in c('metadata','content')){
    report[[elem]] <- list()
  }
  if(!is.null(pcgr_data)){
    report[['metadata']][['pcgr_db_release']] <- pcgr_data[['release_notes']]
    report[['metadata']][['pcgr_version']] <- pcgr_data[['software_version']][['pcgr']]
    report[['metadata']][['cpsr_version']] <- pcgr_data[['software_version']][['cpsr']]
    report[['metadata']][['genome_assembly']] <- pcgr_data[['assembly']][['grch_name']]
    report[['metadata']][['sample_name']] <- sample_name
    report[['metadata']][['report_type']] <- type
    report[['metadata']][['config']] <- config
    report[['metadata']][['tumor_class']] <- 'Not defined'
    report[['metadata']][['tumor_tcga_cohort']] <- 'Not defined'
    report[['metadata']][['tumor_primary_site']] <- 'Not defined'
    report[['metadata']][['medgen_ontology']] <- list()
    report[['metadata']][['medgen_ontology']][['all']] <- pcgr_data[['phenotype_ontology']][['medgen_cancer']]
    report[['metadata']][['medgen_ontology']][['query']] <- NULL

    if(virtual_panel_id >= 0){
      report[['metadata']][['gene_panel']] <- list()
      report[['metadata']][['gene_panel']][['genes']] <- pcgr_data[['virtual_gene_panels']] %>%
        dplyr::filter(id == virtual_panel_id) %>%
        dplyr::select(symbol, confidence_level, panel_name, panel_id, panel_url)
      report[['metadata']][['gene_panel']][['name']] <- unique(report[['metadata']][['gene_panel']][['genes']]$panel_name)
      report[['metadata']][['gene_panel']][['url']] <- unique(report[['metadata']][['gene_panel']][['genes']]$panel_url)
    }

  }

  if(type == 'predisposition'){
    if(!is.null(pcgr_data)){
      report[['metadata']][['medgen_ontology']][['query']] <-
        dplyr::filter(pcgr_data[['phenotype_ontology']][['medgen_cancer']],
                      group == 'Hereditary_Cancer_Susceptibility_NOS' | group == 'Hereditary_Cancer_Syndrome_NOS')
    }
    analysis_element <- 'snv_indel'
    report[['content']][[analysis_element]] <- list()
    report[['content']][[analysis_element]][['eval']] <- FALSE
    report[['content']][[analysis_element]][['variant_display']] <- list()
    report[['content']][[analysis_element]][['variant_set']] <- list()
    report[['content']][[analysis_element]][['zero']] <- FALSE
    for(t in c('class1','class2','class3','class4','class5','gwas','sf')){
      report[['content']][[analysis_element]][['variant_display']][[t]] <- data.frame()
      report[['content']][[analysis_element]][['variant_set']][[t]] <- data.frame()
    }
    report[['content']][[analysis_element]][['variant_set']][['tsv']] <- data.frame()

    for(cl in c('variant_statistic','variant_statistic_cpg','variant_statistic_sf')){
      report[['content']][[analysis_element]][[cl]] <- list()
      for(t in c('n','n_snv','n_indel','n_coding','n_noncoding')){
        report[['content']][[analysis_element]][[cl]][[t]] <- 0
      }
    }

    if(!is.null(report[['metadata']][['config']][['popgen']])){
      if(report[['metadata']][['config']][['popgen']][['pop_gnomad']] != ""){
        pop_tag_info <- pcgrr::get_population_tag(config[['popgen']][['pop_gnomad']], db = "GNOMAD", subset = "non_cancer")
        report[['metadata']][['config']][['popgen']][['vcftag_gnomad']] <- pop_tag_info$vcf_tag
        report[['metadata']][['config']][['popgen']][['popdesc_gnomad']] <- pop_tag_info$pop_description
      }
    }

  }else{

    if(!is.null(pcgr_data)){
      if(config$tumor_type$type != ""){
        report[['metadata']][['tumor_class']] <- config$tumor_type$type
        tumor_group_entry <- dplyr::filter(pcgr_data[['phenotype_ontology']][['cancer_groups']], group == config$tumor_type$type)
        if(nrow(tumor_group_entry) == 1){
          report[['metadata']][['tumor_primary_site']] <- tumor_group_entry$primary_site
          report[['metadata']][['tumor_tcga_cohort']] <- tumor_group_entry$tcga_cohort
          report[['metadata']][['medgen_ontology']][['query']] <- dplyr::filter(pcgr_data[['phenotype_ontology']][['medgen_cancer']], group == config$tumor_type$type)
        }
      }
    }
    for(analysis_element in c('snv_indel','tmb','msi','cna','cna_plot','m_signature',
                              'sequencing_mode','tumor_only','value_box','cna_plot',
                              'tumor_purity','tumor_ploidy','report_display_config')){
      report[['content']][[analysis_element]] <- list()
      report[['content']][[analysis_element]][['eval']] <- FALSE

      if(analysis_element == 'cna_plot'){
        report[['content']][[analysis_element]][['png']] <- NULL
      }
      if(analysis_element == 'report_display_config'){
        report[['content']][[analysis_element]][['opentargets_rank']] <- list()
        report[['content']][[analysis_element]][['opentargets_rank']][['breaks']] <- c(0.40, 0.55, 0.70, 0.85)
        report[['content']][[analysis_element]][['opentargets_rank']][['colors']] <- c("#b8b8ba","#BDD7E7","#6BAED6","#3182BD","#08519C")
      }
      if(analysis_element == 'tumor_purity' | analysis_element == 'tumor_ploidy'){
        report[['content']][[analysis_element]][['estimate']] <- "NA"
        if(!is.null(report[['metadata']][['config']][['tumor_properties']])){
          if(!is.null(report[['metadata']][['config']][['tumor_properties']][[analysis_element]])){
            report[['content']][[analysis_element]][['estimate']] <- report[['metadata']][['config']][['tumor_properties']][[analysis_element]]
            report[['content']][[analysis_element]][['eval']] <- TRUE
          }
        }
      }
      if(analysis_element == 'sequencing_mode'){
        report[['content']][[analysis_element]][['tumor_only']] <- FALSE
        report[['content']][[analysis_element]][['mode']] <- 'Tumor vs. control'
      }

      if(analysis_element == 'value_box'){
        report[['content']][[analysis_element]][['tmb_tertile']] <- 'TMB:\nNot determined'
        report[['content']][[analysis_element]][['msi']] <- 'MSI:\nNot determined'
        report[['content']][[analysis_element]][['scna']] <- 'SCNA:\nNot determined'
        report[['content']][[analysis_element]][['tier1']] <- 'Tier 1 variants:\nNot determined'
        report[['content']][[analysis_element]][['tier2']] <- 'Tier 2 variants:\nNot determined'
        report[['content']][[analysis_element]][['signatures']] <- 'Mutational signatures:\nNot determined'
        report[['content']][[analysis_element]][['tumor_ploidy']] <- 'Tumor ploidy:\nNot provided/determined'
        report[['content']][[analysis_element]][['tumor_purity']] <- 'Tumor purity:\nNot provided/determined'
        report[['content']][[analysis_element]][['tumor_n']] <- 'NA'
      }

      if(analysis_element == 'snv_indel' | analysis_element == 'cna'){
        report[['content']][[analysis_element]][['clinical_evidence_item']] <- list()
        report[['content']][[analysis_element]][['variant_display']] <- list()
        report[['content']][[analysis_element]][['variant_set']] <- list()
        report[['content']][[analysis_element]][['variant_statistic']] <- list()
        report[['content']][[analysis_element]][['variant_statistic']] <- list()
        report[['content']][[analysis_element]][['zero']] <- FALSE
        for(tumorclass in c('any_tumortype','other_tumortype','specific_tumortype')){
          report[['content']][[analysis_element]][['clinical_evidence_item']][[tumorclass]] <- list()
          for(evidence_type in c('prognostic','diagnostic','predictive')){
            for(evidence_level in c('A_B','C_D_E','any')){
              report[['content']][[analysis_element]][['clinical_evidence_item']][[tumorclass]][[evidence_type]][[evidence_level]] <- data.frame()
            }
          }
        }
        if(analysis_element == 'snv_indel'){
          for(t in c('tier1','tier2','tier3','tier4','noncoding')){
            report[['content']][[analysis_element]][['variant_display']][[t]] <- data.frame()
            # if(t == 'tier2' & config$tier_model$tier_model == 'pcgr'){
            #   report[['content']][[analysis_element]][['variant_display']][[t]] <- list()
            #   for(c in c('hotspot','curated_mutation','predicted_driver')){
            #     report[['content']][[analysis_element]][['variant_display']][[t]][[c]] <- data.frame()
            #   }
            # }
            if(t == 'tier3' & config$tier_model$tier_model == 'pcgr_acmg'){
              report[['content']][[analysis_element]][['variant_display']][[t]] <- list()
              for(c in c('proto_oncogene','tumor_suppressor')){
                report[['content']][[analysis_element]][['variant_display']][[t]][[c]] <- data.frame()
              }
            }
          }
          for(t in c('tier1','tier2','tier3','tier4','noncoding','tsv','tsv_unfiltered','maf','coding','all')){
            report[['content']][[analysis_element]][['variant_set']][[t]] <- data.frame()
          }
          for(t in c('n','n_snv','n_indel','n_coding','n_noncoding','n_tier1','n_tier2','n_tier3','n_tier4')){
            report[['content']][[analysis_element]][['variant_statistic']][[t]] <- 0
          }
        }

        if(analysis_element == 'cna'){
          report[['content']][[analysis_element]][['variant_set']][['cna_print']] <- data.frame()
          for(t in c('n_cna_loss','n_cna_gain')){
            report[['content']][[analysis_element]][['variant_statistic']][[t]] <- data.frame()
          }
          for(t in c('segment','oncogene_gain','tsgene_loss','biomarker')){
            report[['content']][[analysis_element]][['variant_display']][[t]] <- data.frame()
          }
          if(config$tier_model$tier_model == 'pcgr_acmg'){
            report[['content']][[analysis_element]][['variant_display']][['biomarkers_tier1']] <- data.frame()
            report[['content']][[analysis_element]][['variant_display']][['biomarkers_tier2']] <- data.frame()
          }
        }
      }
      if(analysis_element == 'm_signature'){
        report[['content']][[analysis_element]][['variant_set']] <- list()
        report[['content']][[analysis_element]][['variant_set']][['all']] <- data.frame()
        report[['content']][[analysis_element]][['missing_data']] <- FALSE
        report[['content']][[analysis_element]][['result']] <- list()
      }
      if(analysis_element == 'tmb'){
        report[['content']][[analysis_element]][['variant_statistic']] <- list()
        report[['content']][[analysis_element]][['variant_statistic']][['n_tmb']] <- 0
        report[['content']][[analysis_element]][['variant_statistic']][['tmb_estimate']] <- 0
        report[['content']][[analysis_element]][['variant_statistic']][['target_size_mb']] <- config[['mutational_burden']][['target_size_mb']]
        report[['content']][[analysis_element]][['variant_statistic']][['tmb_tertile']] <- 'TMB - not determined'
        if(!is.null(pcgr_data)){
          report[['content']][[analysis_element]][['tcga_tmb']] <- pcgr_data[['tcga']][['tmb']]
        }
      }
      if(analysis_element == 'msi'){
        report[['content']][[analysis_element]][['missing_data']] <- FALSE
        report[['content']][[analysis_element]][['prediction']] <- list()
      }
      if(analysis_element == 'tumor_only'){
        report[['content']][[analysis_element]][['variant_set']] <- list()
        report[['content']][[analysis_element]][['variant_set']][['unfiltered']] <- data.frame()
        report[['content']][[analysis_element]][['variant_set']][['filtered']] <- data.frame()
        report[['content']][[analysis_element]][['upset_data']] <- data.frame()
        report[['content']][[analysis_element]][['variant_statistic']] <- list()

        for(successive_filter in c('unfiltered_n','onekg_n_remain','gnomad_n_remain','clinvar_n_remain',
                                   'pon_n_remain','hom_n_remain','het_n_remain','dbsnp_n_remain',
                                   'nonexonic_n_remain','onekg_frac_remain','gnomad_frac_remain','clinvar_frac_remain',
                                   'dbsnp_frac_remain','pon_frac_remain','hom_frac_remain','het_frac_remain',
                                   'nonexonic_frac_remain')){
          report[['content']][[analysis_element]][['variant_statistic']][[successive_filter]] <- 0
        }
      }
    }
  }
  if(!is.null(class)){
    if(!is.null(report[['content']][[class]])){
      return(report[['content']][[class]])
    }
  }

  return(report)
}

#' Function that initates PCGR report object
#'
#' @param report PCGR final report
#' @param report_data Object with PCGR report data
#' @param analysis_element section of PCGR report

update_pcg_report <- function(report, report_data, analysis_element = 'snv_indel'){

  if(!is.null(report_data) & !is.null(report[['content']][[analysis_element]])){
    for(report_elem in names(report_data)){
      if(!is.null(report_data[[report_elem]]) & !is.null(report[['content']][[analysis_element]][[report_elem]])){
         report[['content']][[analysis_element]][[report_elem]] <- report_data[[report_elem]]
      }
    }
  }
  return(report)
}
