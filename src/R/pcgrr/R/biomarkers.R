#' Function that retrieves clinical evidence items (CIVIC, CBMDB) for somatic cancer variants
#'
#' @param sample_calls data frame with sample variants
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_config Object with PCGR configuration parameters
#' @param tumor_type_specificity run against specific tumor type ('tumor_specific', as specified in pcgr_config) or any tumor type 'any'
#' @param biomarker_mapping_stringency - one of 'exact' (allele-specific) or 'approximate' (codon or exon-level biomarkers)
#'
#' @return list
#'
get_clinical_associations_snv_indel <- function(sample_calls, pcgr_data, pcgr_config, tumor_type_specificity = 'any_tumortype', biomarker_mapping_stringency = 1){

  all_eitems <- data.frame()
  variant_set <- data.frame()
  clin_eitems_list <- list()
  for(type in c('diagnostic','prognostic','predictive')){
    clin_eitems_list[[type]] <- list()
    for(elevel in c('any','A_B','C_D_E')){
      clin_eitems_list[[type]][[elevel]] <- data.frame()
    }
  }
  if(tumor_type_specificity == 'any_tumortype'){
    rlogging::message(paste0("Looking up SNV/InDel biomarkers for precision oncology - any tumortype"))
  }else{
    tumor_type_query <- pcgrr::list_to_df(pcgr_config$tumor_type) %>% dplyr::filter(list.element == T) %>% dplyr::select(name)
    if(nrow(tumor_type_query) == 0){
      return(list('clinical_evidence_item' = clin_eitems_list, 'variant_set' = variant_set))
    }
    rlogging::message(paste0("Looking up SNV/InDel biomarkers for precision oncology - ",paste(tumor_type_query$name,collapse=", ")))
  }
  civic_biomarkers <- pcgr_data$civic_biomarkers %>% dplyr::filter(alteration_type == 'MUT')
  cbmdb_biomarkers <- pcgr_data$cbmdb_biomarkers %>% dplyr::filter(alteration_type == 'MUT')
  if("pubmed_html_link" %in% colnames(civic_biomarkers)){
    civic_biomarkers <- dplyr::rename(civic_biomarkers, citation = pubmed_html_link)
  }
  if("evidence_description" %in% colnames(civic_biomarkers)){
    civic_biomarkers <- dplyr::rename(civic_biomarkers, description = evidence_description)
  }
  if("pubmed_html_link" %in% colnames(cbmdb_biomarkers)){
    cbmdb_biomarkers <- dplyr::rename(cbmdb_biomarkers, citation = pubmed_html_link)
  }
  if("evidence_description" %in% colnames(cbmdb_biomarkers)){
    cbmdb_biomarkers <- dplyr::rename(cbmdb_biomarkers, description = evidence_description)
  }
  biomarker_descriptions <- data.frame()

  if(tumor_type_specificity == 'specific_tumortype'){
    tumor_type_query <- pcgrr::list_to_df(pcgr_config$tumor_type) %>% dplyr::filter(list.element == T) %>% dplyr::select(name)
    tumor_tree_query <- dplyr::semi_join(dplyr::select(pcgr_data$phenotype_medgen_oncology,group,do_id,cui,cui_name), tumor_type_query, by=c("group" = "name")) %>% dplyr::filter(!is.na(do_id)) %>% dplyr::distinct()
    civic_biomarkers <- dplyr::semi_join(civic_biomarkers,tumor_tree_query,by=c("disease_ontology_id" = "do_id"))
    cbmdb_biomarkers <- dplyr::semi_join(cbmdb_biomarkers,tumor_tree_query,by=c("disease_ontology_id" = "do_id"))
  }

  #bmarker_mapping_levels <- c('exact','codon')
  bmarker_mapping_levels <- c('exact','codon','exon')
  if(biomarker_mapping_stringency == 2){
    bmarker_mapping_levels <- c('exact','codon','exon','gene')
  }

  for(mapping in bmarker_mapping_levels){
    clinical_evidence_items <- data.frame()
    if(mapping == 'exact'){
      sample_calls_civic <- sample_calls %>% dplyr::filter(!is.na(CIVIC_ID))
      if(nrow(sample_calls_civic) > 0){
        tmp <- dplyr::select(sample_calls_civic,CIVIC_ID,VAR_ID) %>% tidyr::separate_rows(CIVIC_ID,sep=",")
        sample_calls_civic <- merge(tmp,dplyr::select(sample_calls_civic,-c(CIVIC_ID)),by.x = "VAR_ID",by.y = "VAR_ID")
        civic_calls <- dplyr::select(sample_calls_civic,dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
        eitems <- dplyr::inner_join(civic_calls,civic_biomarkers,by=c("CIVIC_ID" = "evidence_id")) %>% dplyr::distinct()
        names(eitems) <- toupper(names(eitems))
        eitems <- eitems %>% dplyr::select(-c(EITEM_CONSEQUENCE,MAPPING_CATEGORY,EITEM_CODON,EITEM_EXON))
        clinical_evidence_items <- rbind(clinical_evidence_items, eitems)
      }
      sample_calls_cbmdb <- sample_calls %>% dplyr::filter(is.na(CIVIC_ID) & !is.na(CBMDB_ID))
      if(nrow(sample_calls_cbmdb) > 0){
        tmp <- dplyr::select(sample_calls_cbmdb,CBMDB_ID,VAR_ID) %>% tidyr::separate_rows(CBMDB_ID,sep=",")
        tmp$CBMDB_ID <- as.integer(tmp$CBMDB_ID)
        sample_calls_cbmdb <- merge(tmp,dplyr::select(sample_calls_cbmdb,-c(CBMDB_ID)),by.x = "VAR_ID",by.y = "VAR_ID")
        cbmdb_calls <- dplyr::select(sample_calls_cbmdb,dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
        eitems <- dplyr::inner_join(cbmdb_calls,cbmdb_biomarkers,by=c("CBMDB_ID" = "evidence_id")) %>% dplyr::distinct()
        names(eitems) <- toupper(names(eitems))
        eitems <- eitems %>% dplyr::select(-c(EITEM_CONSEQUENCE,MAPPING_CATEGORY,EITEM_CODON,EITEM_EXON))
        clinical_evidence_items <- rbind(clinical_evidence_items, eitems)
      }
    }
    else{
      sample_calls_civic <- sample_calls %>% dplyr::filter(!is.na(CIVIC_ID_2))
      if(nrow(sample_calls_civic) > 0){
        tmp <- dplyr::select(sample_calls_civic,CIVIC_ID_2,VAR_ID) %>% tidyr::separate_rows(CIVIC_ID_2,sep=",")
        sample_calls_civic <- merge(tmp,dplyr::select(sample_calls_civic,-c(CIVIC_ID_2)),by.x = "VAR_ID",by.y = "VAR_ID")
        civic_calls <- dplyr::select(sample_calls_civic,dplyr::one_of(pcgr_data$pcgr_all_annotation_columns))
        clinical_evidence_items_all <- dplyr::inner_join(civic_calls,civic_biomarkers,by=c("CIVIC_ID_2" = "evidence_id")) %>% dplyr::distinct()
        names(clinical_evidence_items_all) <- toupper(names(clinical_evidence_items_all))
        if(nrow(clinical_evidence_items_all) > 0){
          if(mapping == 'codon' & 'EITEM_CODON' %in% colnames(clinical_evidence_items_all)){
            if(nrow(clinical_evidence_items_all[!is.na(clinical_evidence_items_all$EITEM_CODON),]) > 0){
              clinical_evidence_items <- dplyr::filter(clinical_evidence_items_all, !is.na(EITEM_CODON))
              clinical_evidence_items <- dplyr::filter(clinical_evidence_items, EITEM_CODON <= AMINO_ACID_END & EITEM_CODON >= AMINO_ACID_START)
              if(nrow(clinical_evidence_items) > 0){
                clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(!is.na(EITEM_CONSEQUENCE) & startsWith(CONSEQUENCE,EITEM_CONSEQUENCE) | is.na(EITEM_CONSEQUENCE))
                if(nrow(clinical_evidence_items) > 0){
                  clinical_evidence_items <- clinical_evidence_items %>% dplyr::select(-c(EITEM_CONSEQUENCE,MAPPING_CATEGORY,EITEM_CODON,EITEM_EXON))
                }else{
                  clinical_evidence_items <- data.frame()
                }
              }else{
                clinical_evidence_items <- data.frame()
              }
            }else{
              clinical_evidence_items <- data.frame()
            }
          }
          if(mapping == 'exon' & 'EITEM_EXON' %in% colnames(clinical_evidence_items_all)){
            if(nrow(clinical_evidence_items_all[!is.na(clinical_evidence_items_all$EITEM_EXON),]) > 0){
              clinical_evidence_items <- dplyr::filter(clinical_evidence_items_all, !is.na(EITEM_EXON))
              clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(EXON == EITEM_EXON)
              if(nrow(clinical_evidence_items) > 0){
                clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(!is.na(EITEM_CONSEQUENCE) & startsWith(CONSEQUENCE,EITEM_CONSEQUENCE) | is.na(EITEM_CONSEQUENCE))
                clinical_evidence_items <- dplyr::filter(clinical_evidence_items, CODING_STATUS == 'coding')
                if(nrow(clinical_evidence_items) > 0){
                  clinical_evidence_items <- clinical_evidence_items %>% dplyr::select(-c(EITEM_CONSEQUENCE,MAPPING_CATEGORY,EITEM_CODON,EITEM_EXON))
                }else{
                  clinical_evidence_items <- data.frame()
                }
              }else{
                clinical_evidence_items <- data.frame()
              }
            }else{
              clinical_evidence_items <- data.frame()
            }
          }
          if(mapping == 'gene' & 'GENESYMBOL' %in% colnames(clinical_evidence_items_all)){
            if(nrow(clinical_evidence_items_all[!is.na(clinical_evidence_items_all$GENESYMBOL),]) > 0){
              clinical_evidence_items <- dplyr::filter(clinical_evidence_items_all, MAPPING_CATEGORY == 'gene')
              if(nrow(clinical_evidence_items) > 0){
                clinical_evidence_items <- clinical_evidence_items %>% dplyr::filter(GENESYMBOL == SYMBOL)
                if(nrow(clinical_evidence_items) > 0){
                  clinical_evidence_items <- dplyr::filter(clinical_evidence_items, CODING_STATUS == 'coding')
                  clinical_evidence_items <- clinical_evidence_items %>% dplyr::select(-c(EITEM_CONSEQUENCE,MAPPING_CATEGORY,EITEM_CODON,EITEM_EXON))
                }else{
                  clinical_evidence_items <- data.frame()
                }
              }else{
                clinical_evidence_items <- data.frame()
              }
            }else{
              clinical_evidence_items <- data.frame()
            }
          }
        }
      }
    }

    if(nrow(clinical_evidence_items) > 0){
      clinical_evidence_items <- clinical_evidence_items %>% dplyr::mutate(EVIDENCE_ID = CIVIC_ID, BIOMARKER_MAPPING = mapping)
      clinical_evidence_items <- dplyr::rename(clinical_evidence_items, HGVSp = HGVSP, HGVSc = HGVSC) %>% dplyr::arrange(EVIDENCE_LEVEL,RATING)
      unique_variants <- clinical_evidence_items %>% dplyr::select(SYMBOL, CONSEQUENCE, CDS_CHANGE, GENOMIC_CHANGE) %>% dplyr::distinct()
      rlogging::message(paste0(nrow(clinical_evidence_items),' clinical evidence item(s) found .. (',nrow(unique_variants),' unique variant(s)), mapping = ',mapping))
      rlogging::message('Underlying variant(s):')
      for(i in 1:nrow(unique_variants)){
        rlogging::message(paste(unique_variants[i,],collapse=" "))
      }
      all_eitems <- rbind(all_eitems, clinical_evidence_items)
    }
    else{
      rlogging::message(paste0(nrow(clinical_evidence_items),' clinical evidence item(s) found .. mapping = ',mapping))
    }

  }

  if(nrow(all_eitems) > 0){
    variant_set <- dplyr::select(all_eitems, dplyr::one_of(pcgr_data$pcgr_all_annotation_columns)) %>%
      dplyr::select(-c(CBMDB_ID,CIVIC_ID,CIVIC_ID_2)) %>% dplyr::distinct()
    clin_eitems_list <- list()
    for(type in c('prognostic','diagnostic','predictive')){
      clin_eitems_list[[type]] <- list()
      clin_eitems_list[[type]][['any']] <- dplyr::select(all_eitems, dplyr::one_of(pcgr_data$tier1_tags_display)) %>% dplyr::filter(EVIDENCE_TYPE == stringr::str_to_title(type)) %>% dplyr::arrange(EVIDENCE_LEVEL)
      if(nrow(clin_eitems_list[[type]][['any']]) > 0){
        if(nrow(clin_eitems_list[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(A|B|B1|B2):"))) > 0){
          clin_eitems_list[[type]][['A_B']] <- clin_eitems_list[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(A|B|B1|B2):")) %>% dplyr::arrange(EVIDENCE_LEVEL)
        }else{
          clin_eitems_list[[type]][['A_B']] <- data.frame()
        }
        if(nrow(clin_eitems_list[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(C|D|E):"))) > 0){
          clin_eitems_list[[type]][['C_D_E']] <- clin_eitems_list[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(C|D|E):")) %>% dplyr::arrange(EVIDENCE_LEVEL)
        }else{
          clin_eitems_list[[type]][['C_D_E']] <- data.frame()
        }
      }else{
        clin_eitems_list[[type]][['C_D_E']] <- data.frame()
        clin_eitems_list[[type]][['A_B']] <- data.frame()
        clin_eitems_list[[type]][['any']] <- data.frame()
      }
    }
  }
  return(list('clinical_evidence_item' = clin_eitems_list, 'variant_set' = variant_set))

}


#' Function that retrieves clinical evidence items (CIVIC, CBMDB) for CNA aberrations
#'
#' @param onco_ts_sets data frame with lost tumor suppressor genes and gained oncogenes
#' @param pcgr_data object with PCGR annotation data
#' @param pcgr_config Object with PCGR configuration parameters
#' @param tumor_type_specificity run against specific tumor type ('tumor_specific', as specified in pcgr_config) or any tumor type 'any'
#'
#' @return list
#'
get_clinical_associations_cna <- function(onco_ts_sets, pcgr_data, pcgr_config, tumor_type_specificity = 'any_tumortype'){

  cna_biomarkers <- data.frame()
  clinical_evidence_items <- list()
  for(type in c('diagnostic','prognostic','predictive')){
    clinical_evidence_items[[type]] <- list()
    for(elevel in c('any','A_B','C_D_E')){
      clinical_evidence_items[[type]][[elevel]] <- data.frame()
    }
  }
  if(tumor_type_specificity == 'any_tumortype'){
    rlogging::message(paste0("Looking up SCNA biomarkers for precision oncology - any tumortype"))
  }else{
    tumor_type_query <- pcgrr::list_to_df(pcgr_config$tumor_type) %>% dplyr::filter(list.element == T) %>% dplyr::select(name)
    if(nrow(tumor_type_query) == 0){
      return(list('clinical_evidence_item' = clinical_evidence_items, 'cna_biomarkers' = cna_biomarkers))
    }
    rlogging::message(paste0("Looking up SCNA biomarkers for precision oncology - ",paste(tumor_type_query$name,collapse=", ")))
  }
  civic_cna_biomarkers <- dplyr::filter(pcgr_data$civic_biomarkers, alteration_type == 'CNA' & !is.na(eitem_consequence)) %>% dplyr::select(genesymbol,evidence_type,evidence_level,evidence_description,cancer_type,evidence_direction,pubmed_html_link,disease_ontology_id,therapeutic_context,rating,clinical_significance,eitem_consequence)
names(civic_cna_biomarkers) <- toupper(names(civic_cna_biomarkers))
  civic_cna_biomarkers <- dplyr::rename(civic_cna_biomarkers, SYMBOL = GENESYMBOL, CNA_TYPE = EITEM_CONSEQUENCE, DESCRIPTION = EVIDENCE_DESCRIPTION, CITATION = PUBMED_HTML_LINK)

  cbmdb_cna_biomarkers <- dplyr::filter(pcgr_data$cbmdb_biomarkers, alteration_type == 'CNA' & !is.na(eitem_consequence)) %>% dplyr::select(genesymbol,evidence_type,evidence_level,evidence_description,cancer_type,evidence_direction,pubmed_html_link,disease_ontology_id, therapeutic_context,rating,clinical_significance,eitem_consequence)
  names(cbmdb_cna_biomarkers) <- toupper(names(cbmdb_cna_biomarkers))
  cbmdb_cna_biomarkers <- dplyr::rename(cbmdb_cna_biomarkers, SYMBOL = GENESYMBOL, CNA_TYPE = EITEM_CONSEQUENCE, DESCRIPTION = EVIDENCE_DESCRIPTION, CITATION = PUBMED_HTML_LINK)

  if(tumor_type_specificity == 'specific_tumortype'){
    tumor_type_query <- pcgrr::list_to_df(pcgr_config$tumor_type) %>% dplyr::filter(list.element == T) %>% dplyr::select(name)
    tumor_tree_query <- dplyr::semi_join(dplyr::select(pcgr_data$phenotype_medgen_oncology,group,do_id,cui,cui_name), tumor_type_query, by=c("group" = "name")) %>% dplyr::filter(!is.na(do_id)) %>% dplyr::distinct()
    civic_cna_biomarkers <- dplyr::semi_join(civic_cna_biomarkers,tumor_tree_query,by=c("DISEASE_ONTOLOGY_ID" = "do_id"))
    cbmdb_cna_biomarkers <- dplyr::semi_join(cbmdb_cna_biomarkers,tumor_tree_query,by=c("DISEASE_ONTOLOGY_ID" = "do_id"))
  }

  for(type in c('tsgene_loss','oncogene_gain')){
    if(!is.null(onco_ts_sets[[type]])){
      if(nrow(onco_ts_sets[[type]]) > 0){
        biomarker_hits <- NULL
        biomarker_hits_civic <- as.data.frame(dplyr::inner_join(onco_ts_sets[[type]], civic_cna_biomarkers, by=c("SYMBOL","CNA_TYPE")))
        if(nrow(biomarker_hits_civic) == 0){
          biomarker_hits_cbmdb <- as.data.frame(dplyr::inner_join(onco_ts_sets[[type]], cbmdb_cna_biomarkers, by=c("SYMBOL","CNA_TYPE")))
          biomarker_hits <- biomarker_hits_cbmdb
        }else{
          biomarker_hits <- biomarker_hits_civic
        }
        if(nrow(biomarker_hits) > 0){
          cna_biomarkers <- rbind(cna_biomarkers,biomarker_hits)
        }
      }
    }
  }

  if(nrow(cna_biomarkers) > 0){
    cna_biomarkers <- cna_biomarkers[c("SYMBOL","CANCER_TYPE","CNA_TYPE","EVIDENCE_LEVEL","CLINICAL_SIGNIFICANCE","EVIDENCE_TYPE","DESCRIPTION","EVIDENCE_DIRECTION","THERAPEUTIC_CONTEXT","CITATION","RATING","GENE_NAME","KEGG_PATHWAY","TARGETED_DRUGS","SEGMENT_LENGTH_MB", "SEGMENT","LogR")]
    cna_biomarkers <- cna_biomarkers %>% dplyr::arrange(EVIDENCE_LEVEL,RATING)
    for(type in c('diagnostic','prognostic','predictive')){
      clinical_evidence_items[[type]] <- list()
      if(nrow(dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == stringr::str_to_title(type))) > 0){
        clinical_evidence_items[[type]][['any']] <- dplyr::filter(cna_biomarkers, EVIDENCE_TYPE == stringr::str_to_title(type))
        if(nrow(clinical_evidence_items[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(A|B|B1|B2):"))) > 0){
          clinical_evidence_items[[type]][['A_B']] <- clinical_evidence_items[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(A|B|B1|B2):"))
        }else{
          clinical_evidence_items[[type]][['A_B']] <- data.frame()
        }
        if(nrow(clinical_evidence_items[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(C|D|E):"))) > 0){
          clinical_evidence_items[[type]][['C_D_E']] <- clinical_evidence_items[[type]][['any']] %>% dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL,"^(C|D|E):"))
        }else{
          clinical_evidence_items[[type]][['C_D_E']] <- data.frame()
        }
      }else{
        for(elevel in c('any','A_B','C_D_E')){
          clinical_evidence_items[[type]][[elevel]] <- data.frame()
        }
      }
    }
  }else{
    for(type in c('diagnostic','prognostic','predictive')){
      clinical_evidence_items[[type]] <- list()
      for(elevel in c('any','A_B','C_D_E')){
        clinical_evidence_items[[type]][[elevel]] <- data.frame()
      }
    }
  }

  return(list('clinical_evidence_item' = clinical_evidence_items, 'cna_biomarkers' = cna_biomarkers))

}
