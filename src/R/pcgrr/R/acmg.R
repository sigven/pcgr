#' Function that tiered variant sets for SNVs/InDels
#'
#' @param sample_calls variant calls subject to mutational signature analysis
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_version PCGR software version
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param genome_seq BSgenome object
#' @param genome_assembly human genome assembly version
#' @param biomarker_mapping_stringency quality level for biomarkers
#' @param callset type of calls
#'
#' @return pcg_report_data data frame with all report elements
#'
generate_report_data_snv_indel_acmg <- function(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config,
                                                genome_seq, genome_assembly ,callset = 'somatic calls', biomarker_mapping_stringency = 1){

  rlogging::message('------')
  rlogging::message(paste0("Generating data for tiered cancer genome report - ",callset, " tier model ",pcgr_config[['tier_model']][['tier_model']],"'"))

  pcg_report_snv_indel <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version , genome_assembly, class = 'snv_indel')

  #coding_consequence_pattern <- "^(stop_|start_lost|frameshift_|missense_variant|splice_donor|splice_acceptor|inframe_)"
  pcg_report_snv_indel[['eval']] <- TRUE
  pcg_report_snv_indel[['variant_set']][['all']] <- sample_calls
  pcg_report_snv_indel[['variant_statistic']][['n']] <- sample_calls %>% nrow()
  pcg_report_snv_indel[['variant_statistic']][['n_snv']] <- sample_calls %>% dplyr::filter(VARIANT_CLASS == 'SNV') %>% nrow()
  pcg_report_snv_indel[['variant_statistic']][['n_indel']] <- sample_calls %>% dplyr::filter(VARIANT_CLASS != 'SNV') %>% nrow()
  pcg_report_snv_indel[['variant_statistic']][['n_coding']] <- sample_calls %>% dplyr::filter(CODING_STATUS == 'coding') %>% nrow()
  rlogging::message(paste0("Number of protein-coding variants: ",pcg_report_snv_indel[['variant_statistic']][['n_coding']]))
  #rlogging::message(paste0("Number of noncoding/silent variants: ",pcg_report_snv_indel[['variant_statistic']][['n_noncoding']]))

  if(pcg_report_snv_indel[['variant_statistic']][['n']] > 0){

    ## Analyze Tier1: Variants of strong clinical significance (Evidence level A+B in tumor type, therapeutic/diagnosis/prognosis)
    ## Analyze Tier2: Variants of potential clinical significance (Evidence level A+B in other tumor types, C+D+E in tumor type, therapeutic/diagnosis/prognosis)
    biomarker_hits_snv_indels_specific <- pcgrr::get_clinical_associations_snv_indel(pcg_report_snv_indel[['variant_set']][['all']],pcgr_data, pcgr_config,tumor_type_specificity = 'specific_tumortype',biomarker_mapping_stringency = biomarker_mapping_stringency)
    biomarker_hits_snv_indels_any <- pcgrr::get_clinical_associations_snv_indel(pcg_report_snv_indel[['variant_set']][['all']],pcgr_data, pcgr_config,tumor_type_specificity = 'any_tumortype',biomarker_mapping_stringency = biomarker_mapping_stringency)
    pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']] <- biomarker_hits_snv_indels_specific$clinical_evidence_item
    pcg_report_snv_indel[['clinical_evidence_item']][['any_tumortype']] <- biomarker_hits_snv_indels_any$clinical_evidence_item
    pcg_report_snv_indel[['variant_set']][['tier1']] <- biomarker_hits_snv_indels_specific$variant_set
    pcg_report_snv_indel[['variant_set']][['tier2']] <- biomarker_hits_snv_indels_any$variant_set


    pcg_report_snv_indel <- pcgrr::assign_tier1_tier2_acmg(pcg_report_snv_indel)
    tier12 <- rbind(pcg_report_snv_indel[['variant_display']][['tier1']],pcg_report_snv_indel[['variant_display']][['tier2']])

    ## Analyze Tier 3: coding mutations in oncogenes/tumor suppressors/cancer census genes
    pcg_report_snv_indel[['variant_set']][['tier3']] <- dplyr::select(pcg_report_snv_indel[['variant_set']][['all']], dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>% dplyr::filter(CODING_STATUS == 'coding') %>% dplyr::filter(ONCOGENE == TRUE | TUMOR_SUPPRESSOR == TRUE)
    if(nrow(tier12) > 0){
      pcg_report_snv_indel[['variant_set']][['tier3']] <- dplyr::anti_join(pcg_report_snv_indel[['variant_set']][['tier3']],tier12, by=c("GENOMIC_CHANGE"))
    }
    tier123 <- tier12
    if(nrow(pcg_report_snv_indel[['variant_set']][['tier3']]) > 0){
      pcg_report_snv_indel[['variant_set']][['tier3']] <- pcg_report_snv_indel[['variant_set']][['tier3']] %>% dplyr::arrange(desc(ONCOSCORE))
      tier123 <- rbind(tier12,dplyr::select(pcg_report_snv_indel[['variant_set']][['tier3']],GENOMIC_CHANGE)) %>% dplyr::distinct()
      pcg_report_snv_indel[['variant_display']][['tier3']][['proto_oncogene']] <- dplyr::select(pcg_report_snv_indel[['variant_set']][['tier3']], dplyr::one_of(pcgr_data$tier2_tags_display)) %>% dplyr::filter(ONCOGENE == TRUE & (is.na(TUMOR_SUPPRESSOR) | TUMOR_SUPPRESSOR == FALSE))
      pcg_report_snv_indel[['variant_display']][['tier3']][['tumor_suppressor']] <- dplyr::select(pcg_report_snv_indel[['variant_set']][['tier3']], dplyr::one_of(pcgr_data$tier2_tags_display)) %>% dplyr::filter(TUMOR_SUPPRESSOR == TRUE)
    }

    ## Analyze Tier 4: Other coding mutations
    pcg_report_snv_indel[['variant_set']][['tier4']] <- dplyr::select(pcg_report_snv_indel[['variant_set']][['all']], dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>% dplyr::filter(CODING_STATUS == 'coding')
    if(nrow(tier123) > 0){
      pcg_report_snv_indel[['variant_set']][['tier4']] <- dplyr::anti_join(pcg_report_snv_indel[['variant_set']][['tier4']],tier123, by=c("GENOMIC_CHANGE"))
    }
    if(nrow(pcg_report_snv_indel[['variant_set']][['tier4']]) > 0){
      pcg_report_snv_indel[['variant_set']][['tier4']] <- pcg_report_snv_indel[['variant_set']][['tier4']] %>% dplyr::arrange(desc(ONCOSCORE))
      pcg_report_snv_indel[['variant_display']][['tier4']] <- dplyr::select(pcg_report_snv_indel[['variant_set']][['tier4']], dplyr::one_of(pcgr_data$tier4_tags_display))
    }

    ## Analyze noncoding mutations
    pcg_report_snv_indel[['variant_set']][['noncoding']] <- dplyr::select(pcg_report_snv_indel[['variant_set']][['all']], dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>% dplyr::filter(CODING_STATUS == 'noncoding')
    if(nrow(pcg_report_snv_indel[['variant_set']][['noncoding']]) > 0){
      pcg_report_snv_indel[['variant_set']][['noncoding']] <- dplyr::anti_join(pcg_report_snv_indel[['variant_set']][['noncoding']],tier123, by=c("GENOMIC_CHANGE"))
      pcg_report_snv_indel[['variant_set']][['noncoding']] <- pcg_report_snv_indel[['variant_set']][['noncoding']] %>% dplyr::arrange(desc(ONCOSCORE))
      pcg_report_snv_indel[['variant_display']][['noncoding']] <- dplyr::select(pcg_report_snv_indel[['variant_set']][['noncoding']], dplyr::one_of(pcgr_data$tier5_tags_display))
    }

    pcg_report_snv_indel[['variant_statistic']][['n_noncoding']] <- pcg_report_snv_indel[['variant_set']][['noncoding']] %>% nrow()
    pcg_report_snv_indel[['variant_set']][['tsv']] <- pcgrr::generate_tier_tsv(pcg_report_snv_indel[['variant_set']], pcgr_data = pcgr_data, pcgr_config, sample_name = sample_name)
    pcg_report_snv_indel[['variant_set']][['maf']] <- pcgrr::tier_to_maf(pcg_report_snv_indel[['variant_set']][['tsv']])
  }

  rlogging::message('------')
  return(pcg_report_snv_indel)

}

#' Function that assigns evidence items for SNVs/InDels to ACMG tiers 1 and 2
#' @param pcg_report_snv_indel report object for snv/indels
#'
#' @return pcg_report_data data frame with all report elements
#'

assign_tier1_tier2_acmg <- function(pcg_report_snv_indel){

  ## get evidence items in other tumor types (A_B) associated with variants that are not found in tier 1
  unique_variants_tier1 <- data.frame()
  unique_variants_tier2 <- data.frame()
  for(type in c('diagnostic','predictive','prognostic')){
    if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']]) > 0){
      vars <- dplyr::select(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']], GENOMIC_CHANGE) %>% dplyr::distinct()
      unique_variants_tier1 <- rbind(unique_variants_tier1, vars) %>% dplyr::distinct()
    }
  }

  for(type in c('diagnostic','predictive','prognostic')){
    if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']]) > 0){
      pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <- pcg_report_snv_indel[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']]
    }
    if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']]) > 0){
      pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <- dplyr::anti_join(pcg_report_snv_indel[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']],pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']],by=c("GENOMIC_CHANGE"))
    }
    if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
      if(nrow(unique_variants_tier1) > 0){
        pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <- dplyr::anti_join(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],unique_variants_tier1,by=c("GENOMIC_CHANGE"))
      }
      if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
        unique_variants_tier2 <- rbind(unique_variants_tier2, dplyr::select(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],GENOMIC_CHANGE)) %>% dplyr::distinct()
      }
    }
    if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']]) > 0){
      if(nrow(unique_variants_tier1) > 0){
        pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']] <- dplyr::anti_join(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']],unique_variants_tier1,by=c("GENOMIC_CHANGE"))
      }
      if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']]) > 0){
        unique_variants_tier2 <- rbind(unique_variants_tier2, dplyr::select(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']],GENOMIC_CHANGE)) %>% dplyr::distinct()
      }
    }
  }


  pcg_report_snv_indel[['variant_display']][['tier1']] <- unique_variants_tier1
  pcg_report_snv_indel[['variant_display']][['tier2']] <- unique_variants_tier2

  if(nrow(unique_variants_tier1) > 0){
    pcg_report_snv_indel[['variant_set']][['tier1']] <- dplyr::semi_join(pcg_report_snv_indel[['variant_set']][['tier1']], unique_variants_tier1, by=c("GENOMIC_CHANGE"))
    pcg_report_snv_indel[['variant_set']][['tier2']] <- dplyr::anti_join(pcg_report_snv_indel[['variant_set']][['tier2']], unique_variants_tier1, by=c("GENOMIC_CHANGE"))
  }
  if(nrow(pcg_report_snv_indel[['variant_set']][['tier2']]) > 0){
    pcg_report_snv_indel[['variant_set']][['tier2']] <- dplyr::semi_join(pcg_report_snv_indel[['variant_set']][['tier2']], unique_variants_tier2, by=c("GENOMIC_CHANGE"))
  }

  return(pcg_report_snv_indel)
}

#' Function that assigns evidence items for SCNAs to ACMG tiers 1 and 2
#' @param pcg_report_snv_indel report object for snv/indels
#'
#' @return pcg_report_data data frame with all report elements
#'


assign_tier1_tier2_acmg_cna <- function(pcg_report_cna){

  ## TIER 1: get evidence items in matching tumor type (A_B level)
  unique_segments_tier1 <- data.frame()
  unique_segments_tier2 <- data.frame()
  for(type in c('diagnostic','predictive','prognostic')){
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']]) > 0){
      vars <- dplyr::select(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']], SEGMENT) %>% dplyr::distinct()
      unique_segments_tier1 <- rbind(unique_segments_tier1, vars) %>% dplyr::distinct()
    }
  }

  ## TIER 2: get evidence items in other tumor types (A_B level) + matching tumor types (C_D_E level) for variants that not found in tier 1
  for(type in c('diagnostic','predictive','prognostic')){
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']]) > 0){
      pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <- pcg_report_cna[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']]
    }
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']]) > 0){
      pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <- dplyr::anti_join(pcg_report_cna[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']],pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']],by=c("SEGMENT"))
    }
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
      if(nrow(unique_segments_tier1) > 0){
        pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <- dplyr::anti_join(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],unique_segments_tier1,by=c("SEGMENT"))
      }
      if(nrow(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
        unique_segments_tier2 <- rbind(unique_segments_tier2, dplyr::select(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],SEGMENT)) %>% dplyr::distinct()
      }
    }
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']]) > 0){
      if(nrow(unique_segments_tier1) > 0){
        pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']] <- dplyr::anti_join(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']],unique_segments_tier1,by=c("SEGMENT"))
      }
      if(nrow(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']]) > 0){
        unique_segments_tier2 <- rbind(unique_segments_tier2, dplyr::select(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']],SEGMENT)) %>% dplyr::distinct()
      }
    }
  }

  pcg_report_cna[['variant_display']][['biomarkers_tier1']] <- unique_segments_tier1
  pcg_report_cna[['variant_display']][['biomarkers_tier2']] <- unique_segments_tier2

  return(pcg_report_cna)
}

#' Function that generates cancer genome report - tier model acmg-like
#'
#' @param project_directory name of project directory
#' @param query_vcf2tsv name of gzipped TSV file (vcf2tsv) with annotated query SNVs/InDels
#' @param pcgr_data List of data frames with PCGR data annotations
#' @param pcg_config Object with PCGR configuration parameters
#' @param cna_segments_tsv name of CNA segments file (tab-separated values)
#' @param sample_name sample identifier
#' @param pcgr_version PCGR software version
#' @param genome_assembly human genome assembly version
#'
#'

generate_report_acmg <- function(project_directory, query_vcf2tsv, pcgr_data, pcgr_config = NULL, sample_name = 'SampleX',
                            cna_segments_tsv = NULL, pcgr_version = '0.6.0', genome_assembly = 'grch37'){

  pcg_report <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = NULL, pcgr_data = pcgr_data)

  genome_seq = BSgenome.Hsapiens.UCSC.hg38
  assembly <- 'hg38'
  if(genome_assembly == 'grch37'){
    genome_seq = BSgenome.Hsapiens.UCSC.hg19
    assembly = 'hg19'
  }

  fnames <- list()
  fnames[['tsv_unfiltered']] <- paste0(project_directory, '/',sample_name,'.pcgr_acmg.snvs_indels.tiers.unfiltered.tsv')
  fnames[['tsv']] <- paste0(project_directory, '/',sample_name,'.pcgr_acmg.snvs_indels.tiers.tsv')
  fnames[['cna_print']] <- paste0(project_directory, '/',sample_name,'.pcgr_acmg.cna_segments.tsv')
  #fnames[['maf']] <- paste0(project_directory, '/',sample_name,'.pcgr.maf')
  fnames[['json']] <- paste0(project_directory, '/',sample_name,'.pcgr_acmg.json')

  if(query_vcf2tsv != 'None'){
    if(!file.exists(query_vcf2tsv) | file.size(query_vcf2tsv) == 0){
      rlogging::warning(paste0("File ",query_vcf2tsv," does not exist or has zero size"))
    }
    else{
      if(!is.null(pcgr_config) & query_vcf2tsv != 'None'){
        sample_calls <- pcgrr::get_calls(query_vcf2tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
        pcg_report_seqmode <- pcgrr::init_pcg_report(pcgr_config, sample_name, pcgr_version, genome_assembly, class = "sequencing_mode")
        pcg_report_seqmode[['eval']] <- TRUE
        if(nrow(sample_calls) > 0){
          if(pcgr_config[['tumor_only']][['vcf_tumor_only']] == TRUE){
            pcg_report_seqmode[['mode']] <- 'Tumor-only (no matching control)'
            pcg_report_seqmode[['tumor_only']] <- TRUE
            pcg_report_tumor_only <- pcgrr::generate_report_data_tumor_only(sample_calls, pcgr_data, pcgr_version,
                                                                            sample_name, pcgr_config, genome_seq, genome_assembly = assembly)
            # pcg_report_snv_indel_unfiltered <- pcgrr::generate_report_data_snv_indel_acmg(pcg_report_tumor_only[['variant_set']][['unfiltered']],
            #                                                                               pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq,
            #                                                                               assembly, callset = "unfiltered callset")
            pcg_report_snv_indel_filtered <- pcgrr::generate_report_data_snv_indel_acmg(pcg_report_tumor_only[['variant_set']][['filtered']],
                                                                                        pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq,
                                                                                        assembly, callset = "germline-filtered callset")

            pcg_report_tumor_only[['variant_set']] <- NULL
            pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_snv_indel_filtered, analysis_element = 'snv_indel')
            pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_tumor_only, analysis_element = 'tumor_only')
          }else{
            pcg_report_snv_indel <- pcgrr::generate_report_data_snv_indel_acmg(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
            pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_snv_indel, analysis_element = 'snv_indel')
            if(pcgr_config[['mutational_signatures']][['mutsignatures']] == T){
              pcg_report_signatures <- pcgrr::generate_report_data_signatures(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
              pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_signatures, analysis_element = 'm_signature')
            }
            if(pcgr_config[['msi']][['msi']] == T){
              pcg_report_msi <- pcgrr::generate_report_data_msi(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
              pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_msi, analysis_element = 'msi')
            }
            if(pcgr_config[['mutational_burden']][['mutational_burden']] == T){
              pcg_report_tmb <- pcgrr::generate_report_data_tmb(sample_calls, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, assembly)
              pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_tmb, analysis_element = 'tmb')
            }
          }
          pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_seqmode, analysis_element = 'sequencing_mode')
        }else{
          pcg_report$snv_indel$zero <- TRUE
          pcg_report[['pcgr_config']][['other']][['list_noncoding']] <- FALSE
          pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_seqmode, analysis_element = 'sequencing_mode')
        }
      }
    }
  }

  if(!is.null(cna_segments_tsv)){
    if(file.exists(cna_segments_tsv)){
      pcg_report_cna <- pcgrr::generate_report_data_cna(cna_segments_tsv, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly, transcript_overlap_pct = pcgr_config[['cna']][['cna_overlap_pct']])
      pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_cna, analysis_element = 'cna')
    }
  }
  for(fname_key in c('tsv','tsv_unfiltered','cna_print')){
    if(fname_key == 'cna_print'){
      if(!is.null(pcg_report[['cna']][['variant_set']][[fname_key]])){
        if(nrow(pcg_report[['cna']][['variant_set']][[fname_key]]) > 0){
          write.table(pcg_report[['cna']][['variant_set']][[fname_key]],file=fnames[[fname_key]], sep="\t",col.names = T,row.names = F,quote = F)
          gzip_command <- paste0('gzip -f ',fnames[[fname_key]])
          system(gzip_command, intern=F)
        }
      }
    }
    else{
      if(!is.null(pcg_report[['snv_indel']][['variant_set']][[fname_key]])){
        if(nrow(pcg_report[['snv_indel']][['variant_set']][[fname_key]]) > 0){
          write.table(pcg_report[['snv_indel']][['variant_set']][[fname_key]],file=fnames[[fname_key]], sep="\t",col.names = T,row.names = F,quote = F)
        }
      }
    }
  }

  pcg_report_value_box <- pcgrr::generate_report_data_value_box(pcg_report, pcgr_data, pcgr_version, sample_name, pcgr_config, genome_seq, genome_assembly)
  pcg_report <- pcgrr::update_pcg_report(pcg_report, pcg_report_value_box, analysis_element = 'value_box')

  for(elem in c('tier1','tier2','tier3','tier4')){
    stat <- paste0('n_',elem)
    pcg_report[['snv_indel']][['variant_statistic']][[stat]] <- nrow(pcg_report[['snv_indel']][['variant_set']][[elem]])
    pcg_report[['snv_indel']][['variant_set']][[elem]] <- NULL
  }
  pcg_report[['snv_indel']][['variant_set']][['noncoding']] <- NULL
  pcg_report[['snv_indel']][['variant_set']][['coding']] <- NULL
  pcg_report[['snv_indel']][['variant_set']][['all']] <- NULL

  pcgr_json <- jsonlite::toJSON(pcg_report, pretty=T,na='string',null = 'null')
  write(pcgr_json, fnames[['json']])
  rmarkdown::render(system.file("templates","report_acmg.Rmd", package="pcgrr"), output_format = rmarkdown::html_document(theme = pcg_report[['pcgr_config']][['visual']][['report_theme']], toc = T, toc_depth = 3, toc_float = T, number_sections = F, includes = rmarkdown::includes(after_body = 'disclaimer.md')), output_file = paste0(sample_name,'.pcgr_acmg.html'), clean = T, output_dir = project_directory, intermediates_dir = project_directory, quiet=T)


}

