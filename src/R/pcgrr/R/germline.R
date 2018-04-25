

#' Function that generates germline-filtered callset and PCGR report statistics for tumor-only callsets
#'
#' @param unfiltered_sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#'
generate_report_data_tumor_only <- function(unfiltered_sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly){

  sample_calls <- unfiltered_sample_calls
  germline_filter_level1_remaining <- 0
  germline_filter_level2_remaining <- 0
  germline_filter_level3_remaining <- 0
  germline_filter_level4_remaining <- 0

  pcg_report_to <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = 'tumor_only')
  sample_calls$dbsnp_germline <- TRUE
  if(nrow(sample_calls[is.na(sample_calls$DBSNPRSID),]) > 0){
    sample_calls[is.na(sample_calls$DBSNPRSID),]$dbsnp_germline <- FALSE
  }
  if(nrow(sample_calls[sample_calls$dbsnp_germline == TRUE & !is.na(sample_calls$DOCM_DISEASE),]) > 0){
    sample_calls[sample_calls$dbsnp_germline == TRUE & !is.na(sample_calls$DOCM_DISEASE),]$dbsnp_germline <- FALSE
  }
  if(nrow(sample_calls[sample_calls$dbsnp_germline == TRUE & (!is.na(sample_calls$CLINVAR_MSID) & !is.na(sample_calls$CLINVAR_VARIANT_ORIGIN) & stringr::str_detect(sample_calls$CLINVAR_VARIANT_ORIGIN,"somatic")),]) > 0){
    sample_calls[sample_calls$dbsnp_germline == TRUE & (!is.na(sample_calls$CLINVAR_MSID) & !is.na(sample_calls$CLINVAR_VARIANT_ORIGIN) & stringr::str_detect(sample_calls$CLINVAR_VARIANT_ORIGIN,"somatic")),]$dbsnp_germline <- FALSE
  }

  sample_calls$tcga_somatic <- TRUE
  if(nrow(sample_calls[is.na(sample_calls$TCGA_FREQUENCY),]) > 0){
    sample_calls[is.na(sample_calls$TCGA_FREQUENCY),]$tcga_somatic <- FALSE
  }
  if(nrow(sample_calls[!is.na(sample_calls$TCGA_PANCANCER_COUNT) & sample_calls$TCGA_PANCANCER_COUNT < pcgr_config[['tumor_only']][['tcga_recurrence']],]) > 0){
    sample_calls[!is.na(sample_calls$TCGA_PANCANCER_COUNT) & sample_calls$TCGA_PANCANCER_COUNT < pcgr_config[['tumor_only']][['tcga_recurrence']],]$tcga_somatic <- FALSE
  }

  rlogging::message(paste0('Total sample calls (', sample_name,'): ',nrow(sample_calls)))
  for(pop in c('EUR','AMR','AFR','SAS','EAS','GLOBAL')){
    sample_calls <- pcgrr::filter_db_germline_variants(sample_calls, pop = pop, dbquery = '1KG', max_tolerated_af = pcgr_config[['tumor_only']][[paste0('maf_onekg_',tolower(pop))]])
  }
  germline_filter_level1_remaining <- nrow(sample_calls)
  rlogging::message(paste0('Excluding coinciding germline variants in 1000 Genomes Project populations'))
  rlogging::message(paste0('Total sample calls remaining: ', nrow(sample_calls)))

  for(pop in c('GLOBAL','NFE','AMR','AFR','SAS','EAS','FIN','OTH')){
    sample_calls <- pcgrr::filter_db_germline_variants(sample_calls, pop = pop, dbquery = 'gnomAD', max_tolerated_af = pcgr_config[['tumor_only']][[paste0('maf_gnomad_',tolower(pop))]])
  }
  germline_filter_level2_remaining <- nrow(sample_calls)
  rlogging::message(paste0('Excluding coinciding germline variants in any population in the genome aggregation database (gnomAD)'))
  rlogging::message(paste0('Total sample calls remaining: ', nrow(sample_calls)))

  germline_filter_level3_remaining <- nrow(sample_calls)
  if(pcgr_config[['tumor_only']][['exclude_dbsnp_nonclinical']] == TRUE){

    if(pcgr_config[['tumor_only']][['keep_known_tcga']] == TRUE){
      rlogging::message(paste0('Excluding non-clinically associated dbSNP variants (dbSNP - not recorded as somatic in DoCM/ClinVar, and not in TCGA with recurrence >= ',pcgr_config[['tumor_only']][['tcga_recurrence']],')'))
      sample_calls <- dplyr::filter(sample_calls, dbsnp_germline == FALSE | (dbsnp_germline == TRUE & tcga_somatic == TRUE))
      rlogging::message(paste0('Total sample calls remaining: ', nrow(sample_calls)))
    }else{
      rlogging::message('Excluding non-clinically associated dbSNP variants (dbSNP - not recorded as somatic in DoCM/ClinVar)')
      sample_calls <- dplyr::filter(sample_calls, dbsnp_germline == FALSE)
      rlogging::message(paste0('Total sample calls remaining: ', nrow(sample_calls)))
    }
    germline_filter_level3_remaining <- nrow(sample_calls)
  }

  germline_filter_level4_remaining <- germline_filter_level3_remaining
  if(pcgr_config[['tumor_only']][['exclude_noncoding']] == TRUE){
    rlogging::message(paste0('Excluding noncoding variants'))
    sample_calls <- dplyr::filter(sample_calls, CODING_STATUS == 'coding')
    rlogging::message(paste0('Total sample calls remaining: ', nrow(sample_calls)))
    germline_filter_level4_remaining <- nrow(sample_calls)
  }

  pcg_report_to[['eval']] <- TRUE
  pcg_report_to[['variant_set']][['unfiltered']] <- unfiltered_sample_calls
  pcg_report_to[['variant_set']][['filtered']] <- sample_calls
  pcg_report_to[['variant_statistic']][['unfiltered_n']] <- nrow(unfiltered_sample_calls)
  pcg_report_to[['variant_statistic']][['onekg_n_remain']] <- germline_filter_level1_remaining
  pcg_report_to[['variant_statistic']][['gnomad_n_remain']] <- germline_filter_level2_remaining
  pcg_report_to[['variant_statistic']][['dbsnp_n_remain']] <- germline_filter_level3_remaining
  pcg_report_to[['variant_statistic']][['noncoding_n_remain']] <- germline_filter_level4_remaining
  for(db_filter in c('onekg','gnomad','dbsnp','noncoding')){
    if(pcg_report_to[['variant_statistic']][[paste0(db_filter,'_n_remain')]] > 0 & pcg_report_to[['variant_statistic']][['unfiltered_n']] > 0){
      pcg_report_to[['variant_statistic']][[paste0(db_filter,'_frac_remain')]] <- round((as.numeric(pcg_report_to[['variant_statistic']][[paste0(db_filter,'_n_remain')]]) /
                                                                                                      pcg_report_to[['variant_statistic']][['unfiltered_n']]) * 100, digits = 2)
    }
  }
  return(pcg_report_to)

}



#' Function that filters a data frame with variants according to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param pop population ('european' or 'global')
#' @param dbquery '1KG' or 'gnomAD'
#' @param min_af minimum allele frequency required for variant to be filtered
#'
#' @return var_df
#'

filter_db_germline_variants <- function(var_df, pop='EUR',dbquery = '1KG', max_tolerated_af = 0.01){

  if(dbquery == '1KG'){
    col <- paste0(pop,"_AF_1KG")
    if(any(grepl(paste0("^",col,"$"),names(var_df)))){
      var_df <- var_df[is.na(var_df[,col]) | var_df[,col] < max_tolerated_af,]
    }
  }
  if(dbquery == 'gnomAD'){
    col <- paste0(pop,"_AF_GNOMAD")
    if(any(grepl(paste0("^",col,"$"),names(var_df)))){
      var_df <- var_df[is.na(var_df[,col]) | var_df[,col] < max_tolerated_af,]
    }
  }

  return(var_df)
}


#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param pop_af_column population_column
#'
#' @return var_df
#'
assign_poplevel_frequency_class <- function(var_df, pop_af_column){

  if(any(grepl(paste0("^",pop_af_column,"$"),names(var_df)))){
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.05,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.05),]$pop_common <- 'Common'
    }
    if(nrow(var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.01 & var_df[pop_af_column] < 0.05),]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.01 & var_df[pop_af_column] < 0.05),]$pop_lowfreq <- 'LowFreq'
    }
    if(nrow(var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.001 & var_df[pop_af_column] < 0.01),]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] >= 0.001 & var_df[pop_af_column] < 0.01),]$pop_rare <- 'Rare'
    }
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] < 0.001 & var_df[pop_af_column] > 0,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] < 0.001 & var_df[pop_af_column] > 0),]$pop_veryrare <- 'VeryRare'
    }
    if(nrow(var_df[!is.na(var_df[pop_af_column]) & var_df[pop_af_column] == 0.00,]) > 0){
      var_df[(!is.na(var_df[pop_af_column]) & var_df[pop_af_column] == 0.00),]$pop_monomorphic <- 'Monomorphic'
    }
  }
  return(var_df)

}

#' Function that assigns a category ('Rare','Common' etc) to population-specific germline frequencies
#'
#' @param var_df data frame with variants
#' @param dbquery 1KG or gnomAD
#' @param pop population
#' @param result_tag name of result column
#'
#' @return var_df
#'
assign_poplevel_frequency <- function(var_df, dbquery='1KG', pop='european', result_tag = 'FREQ_GNOMAD_EUROPEAN'){

  var_df$pop_monomorphic <- rep(NA,nrow(var_df))
  var_df$pop_common <- rep(NA,nrow(var_df))
  var_df$pop_rare <- rep(NA,nrow(var_df))
  var_df$pop_veryrare <- rep(NA,nrow(var_df))
  var_df$pop_lowfreq <- rep(NA,nrow(var_df))

  pop_db <- data.frame('population' = 'american','db' = '1KG', 'tag' = 'AMR_AF_1KG', stringsAsFactors = F)
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = '1KG', 'tag' = 'AFR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'global', 'db' = '1KG', 'tag' = 'GLOBAL_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = '1KG', 'tag' = 'EUR_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = '1KG', 'tag' = 'EAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = '1KG','tag' = 'SAS_AF_1KG'))
  pop_db <- rbind(pop_db, data.frame('population' = 'african', 'db' = 'gnomAD', 'tag' = 'AFR_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'american', 'db' = 'gnomAD', 'tag' = 'AMR_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'european', 'db' = 'gnomAD', 'tag' = 'NFE_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'east_asian', 'db' = 'gnomAD', 'tag' = 'EAS_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'south_asian', 'db' = 'gnomAD','tag' = 'SAS_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'global', 'db' = 'gnomAD','tag' = 'GLOBAL_AF_GNOMAD'))

  tags <- character()
  tags <- c(dplyr::filter(pop_db, population == pop & db == dbquery)$tag)
  for(i in 1:length(tags)){
    var_df <- assign_poplevel_frequency_class(var_df, tags[i])
  }
  var_df[result_tag] <- stringr::str_replace_all(stringr::str_replace_all(paste(var_df$pop_monomorphic,var_df$pop_rare,var_df$pop_veryrare,var_df$pop_lowfreq,var_df$pop_common,sep=","), "(,{0,}NA(,){0,}){1,}",","),"(^,)|(,$)","")
  var_df[nchar(var_df[,result_tag]) == 0,result_tag] <- NA

  var_df <- dplyr::select(var_df, -pop_monomorphic, -pop_common, -pop_rare, -pop_veryrare, -pop_lowfreq)
  return(var_df)
}

