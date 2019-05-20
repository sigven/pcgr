

#' Function that generates germline-filtered callset and PCGR report statistics for a given tumor-only callsets
#'
#' @param unfiltered_sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
generate_report_data_tumor_only <- function(unfiltered_sample_calls, sample_name, pcgr_config){

  sample_calls <- unfiltered_sample_calls
  germline_filter_statistics <- list()
  for(m in c('remain_post_onekg','remain_post_gnomad','remain_post_clinvar','remain_post_dbsnp','remain_post_pon',
             'remain_post_nonexonic','remain_post_hom','remain_post_het')){
    germline_filter_statistics[m] <- 0
  }

  ## initiate report
  pcg_report_to <- pcgrr::init_pcg_report(pcgr_config, sample_name, class = 'tumor_only')

  ## assign evidence tags for germine/somatic state of variants, partially based on user-defined criteria (population allele frequency thresholds)
  sample_calls <- pcgrr::assign_somatic_germline_evidence(sample_calls, pcgr_config)

  ## assign somatic classification based on accumulation of evidence tags and user-defined options
  sample_calls <- pcgrr::assign_somatic_classification(sample_calls, pcgr_config)

  ## Assign statistics to successive filtering levels for different evidence criteria
  ## excluded germline calls found in 1000 Genomes Project
  germline_filter_statistics[['remain_post_onekg']] <- nrow(sample_calls) - nrow(sample_calls[sample_calls$SOMATIC_CLASSIFICATION == "GERMLINE_1KG",])
  rlogging::message(paste0('Excluding coinciding germline variants in 1000 Genomes Project populations'))
  rlogging::message(paste0('Total sample calls remaining: ', germline_filter_statistics[['remain_post_onekg']]))

  ## excluded germline calls found in gnomAD
  germline_filter_statistics[['remain_post_gnomad']] <- germline_filter_statistics[['remain_post_onekg']] - nrow(sample_calls[sample_calls$SOMATIC_CLASSIFICATION == "GERMLINE_GNOMAD",])
  rlogging::message(paste0('Excluding coinciding germline variants in any population in the genome aggregation database (gnomAD)'))
  rlogging::message(paste0('Total sample calls remaining: ', germline_filter_statistics[['remain_post_gnomad']]))

  ## excluded germline calls found in ClinVar
  germline_filter_statistics[['remain_post_clinvar']] <- germline_filter_statistics[['remain_post_gnomad']] - nrow(sample_calls[sample_calls$SOMATIC_CLASSIFICATION == "GERMLINE_CLINVAR",])
  rlogging::message(paste0('Excluding coinciding germline variants in ClinVar'))
  rlogging::message(paste0('Total sample calls remaining: ', germline_filter_statistics[['remain_post_clinvar']]))


  ## excluded germline calls found in panel of normals (if provided)
  germline_filter_statistics[['remain_post_pon']] <- germline_filter_statistics[['remain_post_clinvar']]
  if(pcgr_config[['tumor_only']][['exclude_pon']] == TRUE){
    germline_filter_statistics[['remain_post_pon']] <- germline_filter_statistics[['remain_post_pon']] - nrow(sample_calls[sample_calls$SOMATIC_CLASSIFICATION == "GERMLINE_PON",])
    rlogging::message(paste0('Excluding putative germline variants found in calls from panel-of-normals (PON)'))
    rlogging::message(paste0('Total sample calls remaining: ', germline_filter_statistics[['remain_post_pon']]))
  }

  ## excluded germline calls found with 100% allelic fraction (likely homozygous germline variants)
  germline_filter_statistics[['remain_post_hom']] <- germline_filter_statistics[['remain_post_pon']]
  if(pcgr_config[['tumor_only']][['exclude_likely_hom_germline']] == TRUE){
    germline_filter_statistics[['remain_post_hom']] <- germline_filter_statistics[['remain_post_hom']] - nrow(sample_calls[sample_calls$SOMATIC_CLASSIFICATION == "GERMLINE_HOMOZYGOUS",])
    rlogging::message(paste0('Excluding likely homozygous germline variants found as variants with 100% allelic fraction'))
    rlogging::message(paste0('Total sample calls remaining: ', germline_filter_statistics[['remain_post_hom']]))
  }

  ## excluded germline calls found as likely heterozygous germline variants
  germline_filter_statistics[['remain_post_het']] <- germline_filter_statistics[['remain_post_hom']]
  if(pcgr_config[['tumor_only']][['exclude_likely_het_germline']] == TRUE){
    germline_filter_statistics[['remain_post_het']] <- germline_filter_statistics[['remain_post_het']] - nrow(sample_calls[sample_calls$SOMATIC_CLASSIFICATION == "GERMLINE_HETEROZYGOUS",])
    rlogging::message(paste0('Excluding likely heterozygous germline variants found as variants with 40-60% allelic fraction and recorded in gnomAD + dbSNP'))
    rlogging::message(paste0('Total sample calls remaining: ', germline_filter_statistics[['remain_post_het']]))
  }

  ## excluded calls with dbSNP germline status (if set in config)
  germline_filter_statistics[['remain_post_dbsnp']] <- germline_filter_statistics[['remain_post_het']]
  if(pcgr_config[['tumor_only']][['exclude_dbsnp_nonsomatic']] == TRUE){

    rlogging::message(paste0('Excluding non-somatically associated dbSNP variants (dbSNP - not recorded as somatic in DoCM/ClinVar,
                               and not registered in COSMIC or found in TCGA'))

    germline_filter_statistics[['remain_post_dbsnp']] <- germline_filter_statistics[['remain_post_dbsnp']] - nrow(sample_calls[sample_calls$SOMATIC_CLASSIFICATION == "GERMLINE_DBSNP",])
    rlogging::message(paste0('Total sample calls remaining: ', germline_filter_statistics[['remain_post_dbsnp']]))
  }

  unfiltered_sample_calls <- sample_calls
  sample_calls <- sample_calls %>% dplyr::filter(SOMATIC_CLASSIFICATION == "SOMATIC")

  germline_filter_statistics[['remain_post_nonexonic']] <- germline_filter_statistics[['remain_post_dbsnp']]
  if(pcgr_config[['tumor_only']][['exclude_nonexonic']] == TRUE){
    rlogging::message(paste0('Excluding non-exonic variants'))
    sample_calls <- dplyr::filter(sample_calls, EXONIC_STATUS == "exonic")
    rlogging::message(paste0('Total sample calls remaining: ', nrow(sample_calls)))
    germline_filter_statistics[['remain_post_nonexonic']] <- nrow(sample_calls)
  }


  pcg_report_to[['eval']] <- TRUE
  pcg_report_to[['variant_set']][['unfiltered']] <- unfiltered_sample_calls
  pcg_report_to[['variant_set']][['filtered']] <- sample_calls
  pcg_report_to[['variant_statistic']][['unfiltered_n']] <- nrow(unfiltered_sample_calls)
  pcg_report_to[['variant_statistic']][['onekg_n_remain']] <- germline_filter_statistics[['remain_post_onekg']]
  pcg_report_to[['variant_statistic']][['gnomad_n_remain']] <- germline_filter_statistics[['remain_post_gnomad']]
  pcg_report_to[['variant_statistic']][['clinvar_n_remain']] <- germline_filter_statistics[['remain_post_clinvar']]
  pcg_report_to[['variant_statistic']][['pon_n_remain']] <- germline_filter_statistics[['remain_post_pon']]
  pcg_report_to[['variant_statistic']][['hom_n_remain']] <- germline_filter_statistics[['remain_post_hom']]
  pcg_report_to[['variant_statistic']][['het_n_remain']] <- germline_filter_statistics[['remain_post_het']]
  pcg_report_to[['variant_statistic']][['dbsnp_n_remain']] <- germline_filter_statistics[['remain_post_dbsnp']]
  pcg_report_to[['variant_statistic']][['nonexonic_n_remain']] <- germline_filter_statistics[['remain_post_nonexonic']]
  for(db_filter in c('onekg','gnomad','dbsnp','pon','clinvar','hom','het','nonexonic')){
    if(pcg_report_to[['variant_statistic']][[paste0(db_filter,'_n_remain')]] > 0 & pcg_report_to[['variant_statistic']][['unfiltered_n']] > 0){
      pcg_report_to[['variant_statistic']][[paste0(db_filter,'_frac_remain')]] <- round((as.numeric(pcg_report_to[['variant_statistic']][[paste0(db_filter,'_n_remain')]]) /
                                                                                                      pcg_report_to[['variant_statistic']][['unfiltered_n']]) * 100, digits = 2)
    }
  }
  return(pcg_report_to)

}

#' Function that assigns a SOMATIC_CLASSIFICATION to variants based on evidence found in variant set, potentially limited by user-defined options
#'
#' @param sample_calls data frame with variants
#' @param config configuration object
#'
#' @return sample_calls
#'

assign_somatic_classification <- function(sample_calls, config){

  sample_calls$SOMATIC_CLASSIFICATION <- 'SOMATIC'

  ## Assign non-somatic classification based on various evidence criteria
  ## 1) Frequency of variant in any of the five 1000 Genomes superpopulations is greater than the defined thresholds by the user
  ## 2) Frequency of variant in any of the gnomAD populations is greater than the defined thresholds by the user
  ## 3) Variant is recorded in ClinVar as germline
  ## 4) Variant is found in the user-defined panel-of-normals VCF
  ## 5) Evidence for a likely homozygous germline variant - allelic fraction in tumor sample (AF_TUMOR) is 100% (vary rare for true somatic variants)
  ## 6) Evidence for a likely heterozygous germline variant must satisfy three criteria:
  ##    i) Allelic fraction of alternative allele in tumor sample (AF_TUMOR) is 40-60%,
  ##    ii) Variant is present in dbSNP AND gnomAD
  ##    iii) Variant is neither in COSMIC nor TCGA
  ## 7) Variant is recorded in dbSNP (non-somatic ClinVar/DoCM/COSMIC/TCGA)

  sample_calls <- sample_calls %>%
    dplyr::mutate(SOMATIC_CLASSIFICATION = dplyr::if_else(STATUS_POPFREQ_1KG_ABOVE_TOLERATED == TRUE,"GERMLINE_1KG",SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(SOMATIC_CLASSIFICATION = dplyr::if_else(STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED == TRUE & SOMATIC_CLASSIFICATION == "SOMATIC",
                                                          "GERMLINE_GNOMAD",SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(SOMATIC_CLASSIFICATION = dplyr::if_else(STATUS_CLINVAR_GERMLINE == TRUE & SOMATIC_CLASSIFICATION == "SOMATIC",
                                                          "GERMLINE_CLINVAR",SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(SOMATIC_CLASSIFICATION = dplyr::if_else(STATUS_PON == TRUE & config[['tumor_only']][['exclude_pon']] == TRUE & SOMATIC_CLASSIFICATION == "SOMATIC",
                                                          "GERMLINE_PON",SOMATIC_CLASSIFICATION)) %>%

    dplyr::mutate(SOMATIC_CLASSIFICATION = dplyr::if_else(STATUS_LIKELY_GERMLINE_HOMOZYGOUS == TRUE & config[['tumor_only']][['exclude_likely_hom_germline']] == TRUE
                                                          & SOMATIC_CLASSIFICATION == "SOMATIC","GERMLINE_HOMOZYGOUS",SOMATIC_CLASSIFICATION)) %>%
    dplyr::mutate(SOMATIC_CLASSIFICATION = dplyr::if_else(STATUS_LIKELY_GERMLINE_HETEROZYGOUS == TRUE & config[['tumor_only']][['exclude_likely_het_germline']] == TRUE
                                                          & SOMATIC_CLASSIFICATION == "SOMATIC","GERMLINE_HETEROZYGOUS",SOMATIC_CLASSIFICATION))

  ## set variants found in DBSNP as germline if this option is set to TRUE
  if(config[['tumor_only']][['exclude_dbsnp_nonsomatic']] == TRUE){

    sample_calls <- sample_calls %>%
      dplyr::mutate(SOMATIC_CLASSIFICATION = dplyr::if_else(STATUS_DBSNP_GERMLINE == TRUE & STATUS_TCGA_SOMATIC == FALSE & STATUS_COSMIC == FALSE & SOMATIC_CLASSIFICATION == "SOMATIC","GERMLINE_DBSNP",SOMATIC_CLASSIFICATION))

  }

  return(sample_calls)
}



#' Function that appends several tags denoting evidence for somatic/germline status of variants
#'
#' @param sample_calls data frame with variants
#' @param config configuration object
#'
#' @return sample_calls
#'

assign_somatic_germline_evidence <- function(sample_calls, config){

  ## assign STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED
  for(pop in c('EUR','AMR','AFR','SAS','EAS','GLOBAL')){
    sample_calls <- pcgrr::assign_germline_popfreq_status(sample_calls, pop = pop, dbquery = '1KG', max_tolerated_af = config[['tumor_only']][[paste0('maf_onekg_',tolower(pop))]])
  }

  ## assign STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED
  for(pop in c('GLOBAL','NFE','AMR','AFR','SAS','EAS','ASJ','FIN','OTH')){
    sample_calls <- pcgrr::assign_germline_popfreq_status(sample_calls, pop = pop, dbquery = 'gnomAD', max_tolerated_af = config[['tumor_only']][[paste0('maf_gnomad_',tolower(pop))]])
  }

  ## assign MAX_AF_1KG / MAX_AF_GNOMAD
  gnomad_cols <- c("GLOBAL_AF_GNOMAD", "NFE_AF_GNOMAD", "AMR_AF_GNOMAD", "AFR_AF_GNOMAD", "SAS_AF_GNOMAD",
                   "EAS_AF_GNOMAD", "ASJ_AF_GNOMAD", "FIN_AF_GNOMAD", "OTH_AF_GNOMAD")
  onekg_cols <- c("GLOBAL_AF_1KG", "AMR_AF_1KG", "AFR_AF_1KG", "EAS_AF_1KG", "SAS_AF_1KG", "EUR_AF_1KG")
  sample_calls <- sample_calls %>% dplyr::mutate(MAX_AF_1KG = pmax(!!!rlang::syms(onekg_cols), na.rm = T))
  sample_calls <- sample_calls %>% dplyr::mutate(MAX_AF_GNOMAD = pmax(!!!rlang::syms(gnomad_cols), na.rm = T))

  ## assign STATUS_DBSNP_GERMLINE status to all calls recorded in dbSNP (except relevant in a somatic setting, as defined by ClinVar/DoCM)
  if("DBSNPRSID" %in% colnames(sample_calls)){
    sample_calls <- sample_calls %>%
      dplyr::mutate(STATUS_DBSNP_GERMLINE = dplyr::if_else(!is.na(DBSNPRSID),TRUE,FALSE)) %>%
      dplyr::mutate(STATUS_DBSNP_GERMLINE = dplyr::if_else(STATUS_DBSNP_GERMLINE == T & !is.na(DOCM_PMID),FALSE,STATUS_DBSNP_GERMLINE)) %>%
      dplyr::mutate(STATUS_DBSNP_GERMLINE = dplyr::if_else(STATUS_DBSNP_GERMLINE == T & !is.na(CLINVAR_MSID) & stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"somatic"),FALSE,STATUS_DBSNP_GERMLINE))
  }

  ## assign STATUS_CLINVAR_GERMLINE status to all calls recorded in ClinVar with a "germline" variant-of-origin
  if("CLINVAR_MSID" %in% colnames(sample_calls)){
    sample_calls <- sample_calls %>%
      dplyr::mutate(STATUS_CLINVAR_GERMLINE = dplyr::if_else(!is.na(CLINVAR_MSID) & stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"germline") & !stringr::str_detect(CLINVAR_VARIANT_ORIGIN,"somatic"),TRUE,FALSE))
  }

  ## assign STATUS_LIKELY_GERMLINE_HOMOZYGOUS to all calls with 100% allelic fraction of alternative allele
  if("AF_TUMOR" %in% colnames(sample_calls)){
    sample_calls <- sample_calls %>%
      dplyr::mutate(STATUS_LIKELY_GERMLINE_HOMOZYGOUS = dplyr::if_else(!is.na(AF_TUMOR) & AF_TUMOR == 1,TRUE,FALSE))
  }

  ## assign STATUS_TCGA_SOMATIC to calls in TCGA with recurrence level above user-defined threshold
  if("TCGA_PANCANCER_COUNT" %in% colnames(sample_calls)){
    sample_calls <- sample_calls %>%
      dplyr::mutate(STATUS_TCGA_SOMATIC = dplyr::if_else(!is.na(TCGA_PANCANCER_COUNT),TRUE,FALSE))
  }

  ## assign STATUS_COSMIC to all calls with an identifier in COSMIC
  if("COSMIC_MUTATION_ID" %in% colnames(sample_calls)){
    sample_calls <- sample_calls %>%
      dplyr::mutate(STATUS_COSMIC = dplyr::if_else(!is.na(COSMIC_MUTATION_ID),TRUE,FALSE))
  }

  ## assign STATUS_PON to all calls overlapping the user-defined panel-of-normals VCF
  if("PANEL_OF_NORMALS" %in% colnames(sample_calls)){
    sample_calls <- sample_calls %>%
      dplyr::mutate(STATUS_PON = dplyr::if_else(PANEL_OF_NORMALS == TRUE,TRUE,FALSE))
  }

  ## assign STATUS_LIKELY_GERMLINE_HETEROZYGOUS to all calls that have the alternative allele in the [0.40,0.60] AF range, ii) are registered in dbSNP,
  ## iii) in gnomAD (yet below the user-defined thresholds, and iv) not present in COSMIC/TCGA
  if("AF_TUMOR" %in% colnames(sample_calls) & "MAX_AF_GNOMAD" %in% colnames(sample_calls) & "STATUS_COSMIC" %in% colnames(sample_calls) & "STATUS_TCGA_SOMATIC" %in% colnames(sample_calls)){
    sample_calls <- sample_calls %>%
      dplyr::mutate(STATUS_LIKELY_GERMLINE_HETEROZYGOUS = dplyr::if_else(!is.na(MAX_AF_GNOMAD) & STATUS_DBSNP_GERMLINE == TRUE & !is.na(AF_TUMOR) & AF_TUMOR >= 0.40 & AF_TUMOR <= 0.60 & STATUS_TCGA_SOMATIC == FALSE & STATUS_COSMIC == FALSE,TRUE,FALSE))
  }

  return(sample_calls)
}

#' Function that sets STATUS_POPFREQ_1KG_ABOVE_TOLERATED/STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED to TRUE for variants
#' if any population frequency exceeds max_tolerated_af
#'
#' @param sample_calls data frame with variants
#' @param pop population code (1000 Genomes/gnomAD)
#' @param dbquery 1KG or gnomAD
#' @param max_tolerated_af max tolerated germline allele frequency
#'
#' @return sample_calls
#'
assign_germline_popfreq_status <- function(sample_calls, pop='EUR',dbquery = '1KG', max_tolerated_af = 0.01){


  if(dbquery == '1KG'){
    if(!("STATUS_POPFREQ_1KG_ABOVE_TOLERATED" %in% colnames(sample_calls))){
      sample_calls$STATUS_POPFREQ_1KG_ABOVE_TOLERATED <- FALSE
    }
    col <- paste0(pop,"_AF_1KG")
    if(any(grepl(paste0("^",col,"$"),names(sample_calls)))){
      sample_calls <- sample_calls %>%
        dplyr::mutate(STATUS_POPFREQ_1KG_ABOVE_TOLERATED = dplyr::if_else(!is.na(!!!rlang::sym(col)) & as.numeric(!!!rlang::sym(col)) > max_tolerated_af,TRUE,STATUS_POPFREQ_1KG_ABOVE_TOLERATED))
      }
  }
  if(dbquery == 'gnomAD'){
    if(!("STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED" %in% colnames(sample_calls))){
      sample_calls$STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED <- FALSE
    }
    col <- paste0(pop,"_AF_GNOMAD")
    if(any(grepl(paste0("^",col,"$"),names(sample_calls)))){
      sample_calls <- sample_calls %>%
        dplyr::mutate(STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED = dplyr::if_else(!is.na(!!!rlang::sym(col)) & as.numeric(!!!rlang::sym(col)) > max_tolerated_af,TRUE,STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED))
    }
  }

  return(sample_calls)
}



#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param sample_calls data frame with variants
#' @param pop_af_column population_column
#'
#' @return sample_calls
#'
assign_poplevel_frequency_class <- function(sample_calls, pop_af_column){

  if(any(grepl(paste0("^",pop_af_column,"$"),names(sample_calls)))){
    sample_calls <- sample_calls %>%
      dplyr::mutate(pop_common = dplyr::if_else(!is.na(!!rlang::sym(pop_af_column)) & !!rlang::sym(pop_af_column) >= 0.05,"Common",as.character(NA))) %>%
      dplyr::mutate(pop_lowfreq = dplyr::if_else(!is.na(!!rlang::sym(pop_af_column)) & !!rlang::sym(pop_af_column) >= 0.01 & !!rlang::sym(pop_af_column) < 0.05,"LowFreq",as.character(NA))) %>%
      dplyr::mutate(pop_rare = dplyr::if_else(!is.na(!!rlang::sym(pop_af_column)) & !!rlang::sym(pop_af_column) >= 0.001 & !!rlang::sym(pop_af_column) < 0.01,"Rare",as.character(NA))) %>%
      dplyr::mutate(pop_veryrare = dplyr::if_else(!is.na(!!rlang::sym(pop_af_column)) & !!rlang::sym(pop_af_column) > 0 & !!rlang::sym(pop_af_column) < 0.001,"VeryRare",as.character(NA))) %>%
      dplyr::mutate(pop_monomorphic = dplyr::if_else(!is.na(!!rlang::sym(pop_af_column)) & !!rlang::sym(pop_af_column) == 0,"Monomorphic",as.character(NA)))
  }
  return(sample_calls)

}

#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param sample_calls data frame with variants
#' @param dbquery 1KG or gnomAD
#' @param pop population
#' @param result_tag name of result column
#'
#' @return sample_calls
#'
assign_poplevel_frequency <- function(sample_calls, dbquery='1KG', pop='european', result_tag = 'FREQ_GNOMAD_EUROPEAN'){

  sample_calls$pop_monomorphic <- rep(NA,nrow(sample_calls))
  sample_calls$pop_common <- rep(NA,nrow(sample_calls))
  sample_calls$pop_rare <- rep(NA,nrow(sample_calls))
  sample_calls$pop_veryrare <- rep(NA,nrow(sample_calls))
  sample_calls$pop_lowfreq <- rep(NA,nrow(sample_calls))

  pop_db <- data.frame('population' = 'american','db' = '1KG', 'tag' = 'AMR_AF_1KG', stringsAsFactors = F)
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = '1KG', 'tag' = 'AFR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'global', 'db' = '1KG', 'tag' = 'GLOBAL_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = '1KG', 'tag' = 'EUR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = '1KG', 'tag' = 'EAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = '1KG','tag' = 'SAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = 'gnomAD', 'tag' = 'AFR_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'american', 'db' = 'gnomAD', 'tag' = 'AMR_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = 'gnomAD', 'tag' = 'NFE_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'jewish', 'db' = 'gnomAD', 'tag' = 'ASJ_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'other', 'db' = 'gnomAD', 'tag' = 'OTH_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'finnish', 'db' = 'gnomAD', 'tag' = 'FIN_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = 'gnomAD', 'tag' = 'EAS_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = 'gnomAD','tag' = 'SAS_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'global', 'db' = 'gnomAD','tag' = 'GLOBAL_AF_GNOMAD'))

  tags <- character()
  tags <- c(dplyr::filter(pop_db, population == pop & db == dbquery)$tag)
  for(i in 1:length(tags)){
    sample_calls <- assign_poplevel_frequency_class(sample_calls, tags[i])
  }
  sample_calls[result_tag] <- stringr::str_replace_all(stringr::str_replace_all(paste(sample_calls$pop_monomorphic,sample_calls$pop_rare,sample_calls$pop_veryrare,sample_calls$pop_lowfreq,sample_calls$pop_common,sep=","), "(,{0,}NA(,){0,}){1,}",","),"(^,)|(,$)","")
  sample_calls[nchar(sample_calls[,result_tag]) == 0,result_tag] <- NA

  sample_calls <- dplyr::select(sample_calls, -pop_monomorphic, -pop_common, -pop_rare, -pop_veryrare, -pop_lowfreq)
  return(sample_calls)
}

#' Function that retrieves name of VCF INFO tag and population description for gnomad/1000G population
#'
#' @param population_code three-letter code
#' @param db 1KG or GNOMAD
#' @param subset NA or "non_cancer" (for GNOMAD)
#'
#' @return pop_tag_info
#'
get_population_tag <- function(population_code, db = "1KG", subset = NA){
  pop_tag_info <- list('vcf_tag' = paste0(toupper(population_code),"_AF_",db), "pop_description" = NA)
  if(db == "GNOMAD" & subset == 'non_cancer'){
    pop_tag_info <- list('vcf_tag' = paste0("NON_CANCER_AF_",toupper(population_code)), "pop_description" = NA)
  }

  pop_descriptions_1KG <- data.frame('code' = 'afr', 'pop_description' = 'African', stringsAsFactors = F) %>%
    rbind(data.frame('code' = 'amr', 'pop_description' = 'Admixed American',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'eur', 'pop_description' = 'European',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'eas', 'pop_description' = 'East Asian',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'sas', 'pop_description' = 'South Asian',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'global', 'pop_description' = 'global',stringsAsFactors = F))

  pop_descriptions_gnomad <- data.frame('code' = 'afr', 'pop_description' = 'African', stringsAsFactors = F) %>%
    rbind(data.frame('code' = 'amr', 'pop_description' = 'Admixed American',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'nfe', 'pop_description' = 'Non-Finnish European',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'fin', 'pop_description' = 'Finnish',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'oth', 'pop_description' = 'Other',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'asj', 'pop_description' = 'Ashkenazi Jewish',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'eas', 'pop_description' = 'East Asian',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'sas', 'pop_description' = 'South Asian',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'global', 'pop_description' = 'global',stringsAsFactors = F))

  pop_descriptions_gnomad_non_cancer <- data.frame('code' = 'afr', 'pop_description' = 'African non-cancer subset', stringsAsFactors = F) %>%
    rbind(data.frame('code' = 'amr', 'pop_description' = 'Admixed American non-cancer subset',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'nfe', 'pop_description' = 'Non-Finnish European non-cancer subset',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'fin', 'pop_description' = 'Finnish non-cancer subset',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'oth', 'pop_description' = 'Other non-cancer subset',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'asj', 'pop_description' = 'Ashkenazi Jewish non-cancer subset',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'eas', 'pop_description' = 'East Asian non-cancer subset',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'sas', 'pop_description' = 'South Asian non-cancer subset',stringsAsFactors = F)) %>%
    rbind(data.frame('code' = 'global', 'pop_description' = 'Global non-cancer subset',stringsAsFactors = F))

  if(db == "1KG"){
    pop_entry <- dplyr::filter(pop_descriptions_1KG, code == population_code)
    pop_tag_info[['pop_description']] = pop_entry$pop_description
  }
  if(db == "GNOMAD"){
    pop_entry <- dplyr::filter(pop_descriptions_gnomad, code == population_code)
    pop_tag_info[['pop_description']] = pop_entry$pop_description
    if(subset == 'non_cancer'){
      pop_entry <- dplyr::filter(pop_descriptions_gnomad_non_cancer, code == population_code)
      pop_tag_info[['pop_description']] = pop_entry$pop_description
    }

  }
  return(pop_tag_info)
}

#' Function that makes input data for an UpSet plot (filtering/intersection results) for the somatic-germline
#' classification procedure
#'
#' @param calls unfiltered calls (germline + somatic)
#' @param config config
#'
#' @return upset data
#'
make_upset_plot_data <- function(calls, config){

  columns <- c()
  if(config[['tumor_only']][['exclude_pon']] == TRUE){
    columns <- c(columns,'STATUS_PON')
  }
  if(config[['tumor_only']][['exclude_likely_hom_germline']] == TRUE){
    columns <- c(columns,'STATUS_LIKELY_GERMLINE_HOMOZYGOUS')
  }
  if(config[['tumor_only']][['exclude_likely_het_germline']] == TRUE){
    columns <- c(columns,'STATUS_LIKELY_GERMLINE_HETEROZYGOUS')
  }
  if(config[['tumor_only']][['exclude_dbsnp_nonsomatic']] == TRUE){
    columns <- c(columns,'STATUS_DBSNP_GERMLINE')
  }
  df <- dplyr::select(calls, VAR_ID, STATUS_POPFREQ_1KG_ABOVE_TOLERATED,STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED,STATUS_CLINVAR_GERMLINE)
  for(c in columns){
    if(c %in% colnames(calls)){
      df[,c] <- calls[,c]
    }
  }

  for(v in colnames(df)){
    if(v != 'VAR_ID'){
      df[,v] <- as.integer(df[,v])
    }
  }
  df <- dplyr::rename(df, gnomAD = STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED, OneKGP = STATUS_POPFREQ_1KG_ABOVE_TOLERATED, ClinVar = STATUS_CLINVAR_GERMLINE)
  if("STATUS_PON" %in% colnames(df)){
    df <- dplyr::rename(df, Panel_Of_Normals = STATUS_PON)
  }
  if("STATUS_LIKELY_GERMLINE_HOMOZYGOUS" %in% colnames(df)){
    df <- dplyr::rename(df, HomAF = STATUS_LIKELY_GERMLINE_HOMOZYGOUS)
  }
  if("STATUS_LIKELY_GERMLINE_HETEROZYGOUS" %in% colnames(df)){
    df <- dplyr::rename(df, HetAF = STATUS_LIKELY_GERMLINE_HETEROZYGOUS)
  }
  if("STATUS_DBSNP_GERMLINE" %in% colnames(df)){
    df <- dplyr::rename(df, dbSNP = STATUS_DBSNP_GERMLINE)
  }
  return(df)

}
#' Function that makes an upset calls for germline-filtered variants
#' classification procedure
#'
#' @param upset_data unfiltered calls (germline + somatic)
#'
#' @return p
#'
upset_plot_tumor_only <- function(upset_data){

  isets <- c()
  for(c in colnames(upset_data)){
    if(c != 'VAR_ID'){
      isets <- c(isets,c)
    }
  }
  p <- UpSetR::upset(upset_data, sets = isets, sets.bar.color = "#56B4E9", order.by = "freq", nintersects = 20, text.scale = 1.5, point.size = 6, color.pal = "Blues", empty.intersections = "on")
  return(p)

}

