
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
      pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <-
        dplyr::anti_join(pcg_report_snv_indel[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']],
                         pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']],by=c("GENOMIC_CHANGE"))
    }
    if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
      if(nrow(unique_variants_tier1) > 0){
        pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <-
          dplyr::anti_join(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],unique_variants_tier1,by=c("GENOMIC_CHANGE"))
      }
      if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
        unique_variants_tier2 <- rbind(unique_variants_tier2, dplyr::select(pcg_report_snv_indel[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],GENOMIC_CHANGE)) %>% dplyr::distinct()
      }
    }
    if(nrow(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']]) > 0){
      if(nrow(unique_variants_tier1) > 0){
        pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']] <-
          dplyr::anti_join(pcg_report_snv_indel[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']],unique_variants_tier1,by=c("GENOMIC_CHANGE"))
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
  else{
    pcg_report_snv_indel[['variant_set']][['tier1']] <- data.frame()
  }
  if(nrow(unique_variants_tier2) == 0){
    pcg_report_snv_indel[['variant_set']][['tier2']] <- data.frame()
  }else{
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
      pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <-
        pcg_report_cna[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']]
    }
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']]) > 0){
      pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <-
        dplyr::anti_join(pcg_report_cna[['clinical_evidence_item']][['any_tumortype']][[type]][['A_B']],
                         pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['A_B']],by=c("SEGMENT"))
    }
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
      if(nrow(unique_segments_tier1) > 0){
        pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']] <-
          dplyr::anti_join(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],unique_segments_tier1,by=c("SEGMENT"))
      }
      if(nrow(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']]) > 0){
        unique_segments_tier2 <- rbind(unique_segments_tier2, dplyr::select(pcg_report_cna[['clinical_evidence_item']][['other_tumortype']][[type]][['A_B']],SEGMENT)) %>% dplyr::distinct()
      }
    }
    if(nrow(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']]) > 0){
      if(nrow(unique_segments_tier1) > 0){
        pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']] <-
          dplyr::anti_join(pcg_report_cna[['clinical_evidence_item']][['specific_tumortype']][[type]][['C_D_E']],unique_segments_tier1,by=c("SEGMENT"))
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



