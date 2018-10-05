

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

  for(pop in c('GLOBAL','NFE','AMR','AFR','SAS','EAS','ASJ','FIN','OTH')){
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
  pop_db <- rbind(pop_db, data.frame('population' = 'jewish', 'db' = 'gnomAD', 'tag' = 'ASJ_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'other', 'db' = 'gnomAD', 'tag' = 'OTH_AF_GNOMAD'))
  pop_db <- rbind(pop_db, data.frame('population' = 'finnish', 'db' = 'gnomAD', 'tag' = 'FIN_AF_GNOMAD'))
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

#' Function that retrieves name of VCF INFO tag and population description for gnomad/1000G population
#'
#' @param population_code three-letter code
#' @param db 1KG or GNOMAD
#'
#' @return pop_tag_info
#'
get_population_tag <- function(population_code, db = "1KG"){
  pop_tag_info <- list('vcf_tag' = paste0(toupper(population_code),"_AF_",db), "pop_description" = NA)

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

  if(db == "1KG"){
    pop_entry <- dplyr::filter(pop_descriptions_1KG, code == population_code)
    pop_tag_info[['pop_description']] = pop_entry$pop_description
  }
  if(db == "GNOMAD"){
    pop_entry <- dplyr::filter(pop_descriptions_gnomad, code == population_code)
    pop_tag_info[['pop_description']] = pop_entry$pop_description
  }
  return(pop_tag_info)
}

#' Function that counts insilico predictions of variant effects (i.e. damaging/tolerated) from dbNSFP
#'
#' @param sample_calls sample calls with dbnsfp annotations
#'
#' @return sample_calls
#'
get_insilico_prediction_statistics <- function(sample_calls){

  insilico_pathogenicity_pred_algos <- c('SIFT_DBNSFP','PROVEAN_DBNSFP','META_LR_DBNSFP','FATHMM_DBNSFP','MUTATIONTASTER_DBNSFP',
                                         'MUTATIONASSESSOR_DBNSFP','FATHMM_MKL_DBNSF','M_CAP_DBNSFP','SPLICE_SITE_ADA_DBNSFP','SPLICE_SITE_RF_DBNSFP')
  for(v in c('called','damaging','tolerated','splicing_neutral','splicing_affected')){
    sample_calls[,paste0('n_insilico_',v)] <- 0
  }

  for(algo in insilico_pathogenicity_pred_algos){
    if(algo %in% colnames(sample_calls)){
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] != '.','n_insilico_called'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] != '.','n_insilico_called'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & (sample_calls[,algo] == 'damaging' | sample_calls[,algo] == 'possibly_damaging'),'n_insilico_damaging'] <-  sample_calls[!is.na(sample_calls[,algo]) & (sample_calls[,algo] == 'damaging' | sample_calls[,algo] == 'possibly_damaging'),'n_insilico_damaging'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'tolerated','n_insilico_tolerated'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'tolerated','n_insilico_tolerated'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'affect_splicing','n_insilico_splicing_affected'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'affect_splicing','n_insilico_splicing_affected'] + 1
      sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'splicing_neutral','n_insilico_splicing_neutral'] <-  sample_calls[!is.na(sample_calls[,algo]) & sample_calls[,algo] == 'splicing_neutral','n_insilico_splicing_neutral'] + 1
    }
  }
  return(sample_calls)
}

#' Function that assigns variant pathogenicity scores based on ACMG guidelines
#'
#' @param sample_calls sample calls with dbnsfp annotations
#' @param pcgr_config pcgr configuration object
#' @param pcgr_data pcgr data object
#'
#' @return sample_calls
#'
assign_pathogenicity_score <- function(sample_calls, pcgr_config, pcgr_data){

  path_columns <- c("PVS1","PSC1","PS1","PM1", "PM2", "PP3", "BP4", "n_insilico_called",
                    "n_insilico_damaging", "n_insilico_tolerated", "n_insilico_splicing_neutral",
                    "n_insilico_splicing_affected", "codon_prefix", "clinvar_pathogenic_codon",
                    "clinvar_pathogenic", "hotspot_symbol", "hotspot_codon","hotspot_pvalue",
                    "PATHSCORE","PATHDOC","PATHRANK")
  sample_calls <- sample_calls[, !(colnames(sample_calls) %in% path_columns)]

  sample_calls$PVS1 <- FALSE #PVS1 null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion)
                            #in a gene where LOF is a known mechanism of disease (Dominant mode of inheritance)
  sample_calls$PSC1 <- FALSE #null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion)
                            #in a gene where LOF is a known mechanism of disease (Recessive mode of inheritance)
  sample_calls$PS1 <- FALSE #Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
  #sample_calls$PM1 <- FALSE #Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation
  sample_calls$PM2 <- FALSE #Absent from controls (or at extremely low frequency if recessive) (Table 6) in 1000 Genomes Project, or gnomAD
  sample_calls$PM5 <- FALSE #Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before
  sample_calls$PP3 <- FALSE #Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.)
  sample_calls$BP4 <- FALSE #Multiple lines (>1) of in silico evidence of none deleterious effect.

  ## Assign score based on computational evidence for deleterious/benign effect (predictions from dbNSFP)
  ## Default scheme (from default TOML file):
  ## 1) Damaging: Among 8 possible protein variant effect predictions, at least five algorithms must have made a call, with at least 4 predicted as damaging,
  ##    and at most one predicted as tolerated (PP3)
  ##    - at most 1 prediction for a splicing neutral effect
  ## 2) Tolerated: Among 8 possible protein variant effect predictions, at least five algorithms must have made a call, with at least 4 predicted as tolerated,
  ##    and at most one predicted as damaging (BP4)
  ##    - 0 predictions of splice site affected

  dbnsfp_min_majority <- pcgr_config[['dbnsfp']][['min_majority']]
  dbnsfp_max_minority <- pcgr_config[['dbnsfp']][['max_minority']]
  dbnsfp_min_called <- dbnsfp_max_minority + dbnsfp_min_majority

  sample_calls <- pcgrr::get_insilico_prediction_statistics(sample_calls)
  if(nrow(sample_calls[sample_calls$n_insilico_called >= dbnsfp_min_called &
                       sample_calls$n_insilico_damaging >= dbnsfp_min_majority &
                       sample_calls$n_insilico_tolerated <= dbnsfp_max_minority &
                       sample_calls$n_insilico_splicing_neutral <= 1,]) > 0){
    sample_calls[sample_calls$n_insilico_called >= dbnsfp_min_called &
                   sample_calls$n_insilico_damaging >= dbnsfp_min_majority &
                   sample_calls$n_insilico_tolerated <= dbnsfp_max_minority &
                   sample_calls$n_insilico_splicing_neutral <= 1,]$PP3 <- TRUE
  }
  if(nrow(sample_calls[sample_calls$n_insilico_called >= dbnsfp_min_called &
                       sample_calls$n_insilico_tolerated >= dbnsfp_min_majority &
                       sample_calls$n_insilico_damaging <= dbnsfp_max_minority &
                       sample_calls$n_insilico_splicing_affected == 0,]) > 0){
    sample_calls[sample_calls$n_insilico_called >= dbnsfp_min_called &
                   sample_calls$n_insilico_tolerated >= dbnsfp_min_majority &
                   sample_calls$n_insilico_damaging <= dbnsfp_max_minority &
                   sample_calls$n_insilico_splicing_affected == 0,]$BP4 <- TRUE
  }
  if(nrow(sample_calls[sample_calls$n_insilico_splicing_affected == 2,]) > 0){
    sample_calls[sample_calls$n_insilico_splicing_affected == 2,]$PP3 <- TRUE
  }

  ## Assign score based on absence/extremely low frequency from controls (1000G/gnomAD) - PM2
  if('GLOBAL_AF_1KG' %in% colnames(sample_calls) & 'GLOBAL_AF_GNOMAD' %in% colnames(sample_calls)){
    if(nrow(sample_calls[is.na(sample_calls$GLOBAL_AF_1KG) & is.na(sample_calls$GLOBAL_AF_GNOMAD),]) > 0){
      sample_calls[is.na(sample_calls$GLOBAL_AF_1KG) & is.na(sample_calls$GLOBAL_AF_GNOMAD),]$PM2 <- TRUE
    }
    if(nrow(sample_calls[!is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD < 0.0005,]) > 0){
      sample_calls[!is.na(sample_calls$GLOBAL_AF_GNOMAD) & sample_calls$GLOBAL_AF_GNOMAD < 0.0005,]$PM2 <- TRUE
    }
  }
  ## Assign score based on loss-of-function variant in known predisposition gene
  gain_of_function_genes <- c('ALK','SOS1','MET','EGFR','RET','HRAS','CDK4','KIT','PDGFRA','RHBDF2','PTPN11')

  if(nrow(sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI),])){
    if(nrow(sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"Dom|Dom&Rec"),]) > 0){
      sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"Dom|Dom&Rec"),]$PVS1 <- TRUE
    ## skip genes that are associated with gain-of-function
      if(nrow(sample_calls[sample_calls$PVS1 == TRUE & sample_calls$SYMBOL %in% gain_of_function_genes,]) > 0){
        sample_calls[sample_calls$PVS1 == TRUE & sample_calls$SYMBOL %in% gain_of_function_genes,]$PVS1 <- FALSE
      }
    }

    if(nrow(sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"Rec"),]) > 0){
      sample_calls[!is.na(sample_calls$LOSS_OF_FUNCTION) & sample_calls$LOSS_OF_FUNCTION == TRUE & !is.na(sample_calls$CANCER_PREDISPOSITION_MOI) & stringr::str_detect(sample_calls$CANCER_PREDISPOSITION_MOI,"Rec"),]$PSC1 <- TRUE
      ## skip genes that are associated with gain-of-function
      if(nrow(sample_calls[sample_calls$PSC1 == TRUE & sample_calls$SYMBOL %in% gain_of_function_genes,]) > 0){
        sample_calls[sample_calls$PSC1 == TRUE & sample_calls$SYMBOL %in% gain_of_function_genes,]$PSC1 <- FALSE
      }
    }
  }

  ## Assign score (PS1/PM5) to variants that are
  ## 1) coinciding with known pathogenic missense variants (yet with different nucleotide change) - PS1
  ## 2) occurs at the same codon as a known pathogenic missense variant - PM5
  sample_calls$codon_prefix <- NA
  if(nrow(sample_calls[!is.na(sample_calls$CONSEQUENCE) & sample_calls$CONSEQUENCE == 'missense_variant',]) > 0){
    sample_calls[sample_calls$CONSEQUENCE == 'missense_variant',]$codon_prefix <- stringr::str_match(sample_calls[sample_calls$CONSEQUENCE == 'missense_variant',]$HGVSp_short,"p\\.[A-Z]{1}[0-9]{1,}")
  }
  if(nrow(sample_calls[!is.na(sample_calls$codon_prefix),]) > 0){
    sample_calls_pathogenic_codon <- dplyr::left_join(dplyr::filter(dplyr::select(sample_calls,VAR_ID,codon_prefix,SYMBOL), !is.na(codon_prefix)), pcgr_data$clinvar_pathogenic_cpg_codon, by=c("codon_prefix" = "codon_prefix","SYMBOL" = "symbol"))
    sample_calls <- dplyr::left_join(sample_calls, dplyr::select(sample_calls_pathogenic_codon, VAR_ID, clinvar_pathogenic_codon), by=c("VAR_ID"))
  }
  if(nrow(sample_calls[!is.na(sample_calls$HGVSp_short),]) > 0){
    sample_calls_pathogenic_hgvsp <- dplyr::left_join(dplyr::filter(dplyr::select(sample_calls, VAR_ID, HGVSp_short, SYMBOL), !is.na(HGVSp_short)), pcgr_data$clinvar_pathogenic_cpg_hgvsp, by=c("HGVSp_short" = "hgvs_p","SYMBOL" = "symbol"))
    sample_calls <- dplyr::left_join(sample_calls, dplyr::select(sample_calls_pathogenic_hgvsp, VAR_ID, clinvar_pathogenic), by=c("VAR_ID"))
  }

  if(nrow(sample_calls[!is.na(sample_calls$clinvar_pathogenic) & sample_calls$clinvar_pathogenic == T & (stringr::str_detect(sample_calls$CLINVAR_CLNSIG,"uncertain_significance|not_provided") | is.na(sample_calls$CLINVAR_MSID)),]) > 0){
    sample_calls[!is.na(sample_calls$clinvar_pathogenic) & sample_calls$clinvar_pathogenic == T &  (stringr::str_detect(sample_calls$CLINVAR_CLNSIG,"uncertain_significance|not_provided") | is.na(sample_calls$CLINVAR_MSID)),]$PS1 <- TRUE
  }
  if(nrow(sample_calls[!is.na(sample_calls$clinvar_pathogenic_codon) & sample_calls$clinvar_pathogenic_codon == T &  (stringr::str_detect(sample_calls$CLINVAR_CLNSIG,"uncertain_significance|not_provided")| is.na(sample_calls$CLINVAR_MSID)),]) > 0){
    sample_calls[!is.na(sample_calls$clinvar_pathogenic_codon) & sample_calls$clinvar_pathogenic_codon == T & (stringr::str_detect(sample_calls$CLINVAR_CLNSIG,"uncertain_significance|not_provided")| is.na(sample_calls$CLINVAR_MSID)),]$PM5 <- TRUE
  }

  ## if previously found coinciding with pathogenic variant, set PM5 to false
  if(nrow(sample_calls[is.na(sample_calls$PS1) & is.na(sample_calls$PM5) & sample_calls$PM5 == T & sample_calls$PS1 == T,]) > 0){
    sample_calls[is.na(sample_calls$PS1) & is.na(sample_calls$PM5) & sample_calls$PM5 == T & sample_calls$PS1 == T,]$PM5 <- FALSE
  }

  sample_calls <- sample_calls %>% tidyr::separate(CANCER_MUTATION_HOTSPOT, c("hotspot_symbol","hotspot_codon","hotspot_pvalue"),sep="\\|",remove=F,extra="drop")
  if(nrow(sample_calls[!is.na(sample_calls$hotspot_codon),]) > 0){
    sample_calls[!is.na(sample_calls$hotspot_codon),]$hotspot_codon <- paste0('p.',sample_calls[!is.na(sample_calls$hotspot_codon),]$hotspot_codon)
  }

  sample_calls_somatic_hotspots <- sample_calls %>%
    dplyr::filter(!is.na(hotspot_codon) & !is.na(hotspot_symbol)) %>%
    dplyr::filter(!is.na(codon_prefix) & !is.na(SYMBOL)) %>%
    dplyr::filter(SYMBOL == hotspot_symbol & hotspot_codon == codon_prefix) %>%
    dplyr::select(VAR_ID) %>%
    dplyr::distinct() %>%
    dplyr::mutate(PM1 = TRUE)

  if(nrow(sample_calls_somatic_hotspots) > 0){
    sample_calls <- dplyr::left_join(sample_calls, sample_calls_somatic_hotspots, by=c("VAR_ID"))
  }else{
    sample_calls$PM1 <- FALSE
  }
  if(nrow(sample_calls[is.na(sample_calls$PM1),]) > 0){
    sample_calls[is.na(sample_calls$PM1),]$PM1 <- FALSE
  }

  sample_calls$PATHSCORE <- 0
  sample_calls$PATHDOC <- ""
  if(nrow(sample_calls[sample_calls$PVS1 == TRUE,]) > 0){
    sample_calls[sample_calls$PVS1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PVS1 == TRUE,]$PATHSCORE + 8
    sample_calls[sample_calls$PVS1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PVS1 == TRUE,]$PATHDOC,"- Loss-of-function variant in known susceptibility/syndrome gene (dominant MOI)",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PSC1 == TRUE,]) > 0){
    sample_calls[sample_calls$PSC1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PSC1 == TRUE,]$PATHSCORE + 4
    sample_calls[sample_calls$PSC1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PSC1 == TRUE,]$PATHDOC,"- Loss-of-function variant in known susceptibility/syndrome gene (recessive MOI)",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PS1 == TRUE,]) > 0){
    sample_calls[sample_calls$PS1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PS1 == TRUE,]$PATHSCORE + 7
    sample_calls[sample_calls$PS1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PS1 == TRUE,]$PATHDOC,"- Same peptide change as a previously established pathogenic variant",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PM1 == TRUE,]) > 0){
    sample_calls[sample_calls$PM1 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PM1 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PM1 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PM1 == TRUE,]$PATHDOC,"- Variant located in somatic mutation hotspot",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PM2 == TRUE,]) > 0){
    sample_calls[sample_calls$PM2 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PM2 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PM2 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PM2 == TRUE,]$PATHDOC,"- Absent or extremely low frequency in the general population (gnomAD global MAF < 0.0005)",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PM5 == TRUE,]) > 0){
    sample_calls[sample_calls$PM5 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PM5 == TRUE,]$PATHSCORE + 2
    sample_calls[sample_calls$PM5 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PM5 == TRUE,]$PATHDOC,"- Different peptide change of a pathogenic variant at the same reference peptide",sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$PP3 == TRUE,]) > 0){
    sample_calls[sample_calls$PP3 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$PP3 == TRUE,]$PATHSCORE + 1
    sample_calls[sample_calls$PP3 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$PP3 == TRUE,]$PATHDOC,paste0("- Multiple lines (>=",pcgr_config[['dbnsfp']][['min_majority']],") of in silico evidence of deleterious effect"),sep="<br>")
  }
  if(nrow(sample_calls[sample_calls$BP4 == TRUE,]) > 0){
    sample_calls[sample_calls$BP4 == TRUE,]$PATHSCORE <- sample_calls[sample_calls$BP4 == TRUE,]$PATHSCORE - 1
    sample_calls[sample_calls$BP4 == TRUE,]$PATHDOC <- paste(sample_calls[sample_calls$BP4 == TRUE,]$PATHDOC,paste0("- Multiple lines (>=",pcgr_config[['dbnsfp']][['min_majority']],") of in silico evidence of benign effect"),sep="<br>")
  }

  sample_calls$PATHDOC <- stringr::str_replace(sample_calls$PATHDOC,"^<br>","")
  sample_calls$PATHDOC <- paste0("CPSR pathogenicity score: ",sample_calls$PATHSCORE, "<br>",sample_calls$PATHDOC)

  sample_calls$PATHRANK <- NA
  if(nrow(sample_calls[sample_calls$PATHSCORE > 8,]) > 0){
    sample_calls[sample_calls$PATHSCORE > 8,]$PATHRANK <- 'HIGH'
  }
  if(nrow(sample_calls[sample_calls$PATHSCORE > 4 & sample_calls$PATHSCORE <= 8,]) > 0){
    sample_calls[sample_calls$PATHSCORE > 4 & sample_calls$PATHSCORE <= 8,]$PATHRANK <- 'MODERATE'
  }
  if(nrow(sample_calls[sample_calls$PATHSCORE <= 4,]) > 0){
    sample_calls[sample_calls$PATHSCORE <= 4,]$PATHRANK <- 'LOW'
  }
  return(sample_calls)

}



#' Function that generates tiered annotated variant datasets for CPSR
#'
#' @param pcg_report List with tiered variants
#'
#' @return tsv_variants data frame with tier-annotated list of variants for tab-separated output
#'
generate_tier_tsv_cpsr <- function(pcg_report, sample_name = "test"){

  predispose_tsv_tags <- c("VAR_ID","VCF_SAMPLE_ID","CODING_STATUS","SYMBOL","ENSEMBL_GENE_ID","ENSEMBL_TRANSCRIPT_ID","GENOTYPE","CONSEQUENCE","PROTEIN_CHANGE",
                           "GENE_NAME","PROTEIN_DOMAIN", "HGVSp","HGVSc", "CDS_CHANGE","PROTEIN_FEATURE","EFFECT_PREDICTIONS",
                           "LOSS_OF_FUNCTION", "DBSNP","CLINVAR_CLINICAL_SIGNIFICANCE", "CLINVAR_MSID","CLINVAR_VARIANT_ORIGIN",
                           "CLINVAR_CONFLICTED", "CLINVAR_PHENOTYPE","VEP_ALL_CONSEQUENCE", "ONCOGENE", "ONCOSCORE","TUMOR_SUPPRESSOR",
                           "PATHSCORE","PATHRANK", "PATHDOC","GLOBAL_AF_GNOMAD", pcg_report[['pcgr_config']][['popgen']][['vcftag_gnomad']],
                           "GLOBAL_AF_1KG", pcg_report[['pcgr_config']][['popgen']][['vcftag_tgp']],"TIER","TIER_DESCRIPTION",
                           "GENOMIC_CHANGE", "GENOME_VERSION")

  rlogging::message("Generating tiered set of result variants for output in tab-separated values (TSV) file")

  tsv_variants <- data.frame()
  for(tier in c("tier1", "tier2", "tier3A", "tier3B","gwas")){
    if(tier != 'tier3B' & tier != "gwas"){
      predispose_tags <- predispose_tsv_tags
      for(ph in c('cancer_phenotype','noncancer_phenotype')){
        tierset <- data.frame()
        if(nrow(pcg_report[['snv_indel']][['variant_display']][[tier]][[ph]]) > 0){
          tierset <- pcg_report[['snv_indel']][['variant_display']][[tier]][[ph]]
          tierset$VCF_SAMPLE_ID <- sample_name
          if(tier == 'tier1'){
            tierset$TIER <- 'TIER 1'
            if(ph == 'cancer_phenotype'){
              tierset$TIER_DESCRIPTION <- 'Tier 1 - Pathogenic variant (ClinVar) associated with cancer phenotype'
            }else{
              tierset$TIER_DESCRIPTION <- 'Tier 1 - Pathogenic variant (ClinVar) associated with undefined/noncancer phenotype'
            }
          }
          if(tier == 'tier2'){
            tierset$TIER <- 'TIER 2'
            if(ph == 'cancer_phenotype'){
              tierset$TIER_DESCRIPTION <- 'Tier 2 - Likely pathogenic variant (ClinVar) associated with cancer phenotype'
            }else{
              tierset$TIER_DESCRIPTION <- 'Tier 2 - Likely pathogenic variant (ClinVar) associated with undefined/noncancer phenotype'
            }
          }
          if(tier == 'tier3A'){
            tierset$TIER <- 'TIER 3A'
            if(ph == 'cancer_phenotype'){
              tierset$TIER_DESCRIPTION <- 'Tier 3A - Variant of uncertain significance (VUS in ClinVar) associated with cancer phenotype'
            }else{
              tierset$TIER_DESCRIPTION <- 'Tier 3A - Variant of uncertain significance (VUS in ClinVar) associated with undefined/noncancer phenotype'
            }
          }
        }
        if(!is.null(pcg_report[['pcgr_config']][['custom_tags']])){
          if(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']] != ""){
            tags <- stringr::str_split(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']],pattern = ",")[[1]]
            for(t in tags){
              t <- stringr::str_trim(t)
              if(t %in% colnames(vcf_data_df)){
                predispose_tags <- c(predispose_tsv_tags,t)
              }
            }
          }
        }
        if(nrow(tierset) > 0){
          tsv_variants <- as.data.frame(dplyr::bind_rows(tsv_variants, dplyr::select(tierset, dplyr::one_of(predispose_tags))))
        }
      }
    }else{
      predispose_tags <- predispose_tsv_tags
      if(!is.null(pcg_report[['pcgr_config']][['custom_tags']])){
        if(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']] != ""){
          tags <- stringr::str_split(pcg_report[['pcgr_config']][['custom_tags']][['custom_tags']],pattern = ",")[[1]]
          for(t in tags){
            t <- stringr::str_trim(t)
            if(t %in% colnames(vcf_data_df)){
              predispose_tags <- c(predispose_tsv_tags,t)
            }
          }
        }
      }
      if(nrow(pcg_report[['snv_indel']][['variant_display']][[tier]]) > 0){
        tierset <- pcg_report[['snv_indel']][['variant_display']][[tier]]
        tierset$VCF_SAMPLE_ID <- sample_name
        if(tier == 'tier3B'){
          tierset$TIER_DESCRIPTION <- 'Tier 3B - Unclassified variants'
          tierset$TIER <- 'TIER 3B'
        }
        if(tier == 'gwas'){
          tierset$TIER_DESCRIPTION <- 'GWAS hit'
          tierset$TIER <- 'GWAS'
        }
        if(nrow(tierset) > 0){
          tsv_variants <- as.data.frame(dplyr::bind_rows(tsv_variants, dplyr::select(tierset, dplyr::one_of(predispose_tags))))
        }
      }
    }

  }
  tsv_variants$DBSNP <- unlist(lapply(stringr::str_match_all(tsv_variants$DBSNP,">rs[0-9]{1,}<"),paste,collapse=","))
  tsv_variants$DBSNP <- stringr::str_replace_all(tsv_variants$DBSNP,">|<", "")
  tsv_variants$GENE_NAME <- unlist(lapply(stringr::str_match_all(tsv_variants$GENE_NAME,">.+<"),paste,collapse=","))
  tsv_variants$GENE_NAME <- stringr::str_replace_all(tsv_variants$GENE_NAME,">|<", "")
  tsv_variants$PROTEIN_DOMAIN <- unlist(lapply(stringr::str_match_all(tsv_variants$PROTEIN_DOMAIN,">.+<"),paste,collapse=","))
  tsv_variants$PROTEIN_DOMAIN <- stringr::str_replace_all(tsv_variants$PROTEIN_DOMAIN,">|<", "")
  tsv_variants <- tsv_variants %>% dplyr::distinct()

  pcg_report[['snv_indel']][['variant_set']][['tsv']] <- tsv_variants
  for(t in c('TIER 1','TIER 2','TIER 3A','TIER 3B','GWAS')){
    if(t == 'TIER 1'){
      pcg_report[['snv_indel']][['variant_set']][['tier1']] <- dplyr::filter(tsv_variants, TIER == 'TIER 1')
    }
    if(t == 'TIER 2'){
      pcg_report[['snv_indel']][['variant_set']][['tier2']] <- dplyr::filter(tsv_variants, TIER == 'TIER 2')
    }
    if(t == 'TIER 3A'){
      pcg_report[['snv_indel']][['variant_set']][['tier3A']] <- dplyr::filter(tsv_variants, TIER == 'TIER 3A')
    }
    if(t == 'TIER 3B'){
      pcg_report[['snv_indel']][['variant_set']][['tier3B']] <- dplyr::filter(tsv_variants, TIER == 'TIER 3B')
    }
    if(t == 'GWAS'){
      pcg_report[['snv_indel']][['variant_set']][['gwas']] <- dplyr::filter(tsv_variants, TIER == 'GWAS')
    }
  }
  return(pcg_report)
}

