#' Function that retrieves clinical evidence items (CIVIC, CBMDB) for somatic cancer variants
#'
#' @param sample_calls data frame with sample variants
#' @param pcgr_data object with PCGR annotation data
#' @param tumor_type Primary site of tumor sample
#' @param tumor_type_specificity run against specific tumor type ('tumor_specific', as specified in pcgr_config) or any tumor type 'any'
#' @param biomarker_mapping_stringency - one of 'exact' (allele-specific) or 'approximate' (codon or exon-level biomarkers)
#'
#' @return list

get_clin_assocs_snv_indel <- function(sample_calls, pcgr_data, tumor_type,
                                      tumor_type_specificity = "any",
                                      mapping_stringency = 1) {

  invisible(assertthat::assert_that(tumor_type_specificity == "any" | tumor_type_specificity == "specific",
                          msg = paste0("Argument 'tumor_type_specificy' can only take two values: 'any' or 'specific' and NOT: ",
                                       tumor_type_specificity)))
  if(tumor_type_specificity == "specific"){
    invisible(assertthat::assert_that(tumor_type != "Cancer, NOS",
                                    msg = paste0("If tumor_specificy is 'specific', 'tumor_type' must be specified and NOT:  ",
                                                 tumor_type)))
  }

  ## Initialize lists that hold clinical evidence items - across evidence types (diagnostic/prognostic/predictive)
  ## and across evidence levels
  all_var_evidence_items <- data.frame()
  variant_set <- data.frame()
  clin_eitems_list <- list()
  for (type in c("diagnostic", "prognostic", "predictive")) {
    clin_eitems_list[[type]] <- list()
    for (elevel in c("any", "A_B", "C_D_E")) {
      clin_eitems_list[[type]][[elevel]] <- data.frame()
    }
  }

  ## load all somatic clinical evidence items/biomarkers (civic, cgi)
  biomarkers <- pcgrr::load_all_eitems(pcgr_data = pcgr_data, alteration_type = "MUT",
                                       origin = "Somatic")

  ## limit clinical evidence items by primary tumor site
  if (tumor_type_specificity == "any") {
    rlogging::message(paste0("Looking up SNV/InDel biomarkers for precision oncology - any tumortype"))
  }else{

    biomarkers <-
      pcgrr::filter_eitems_by_site(biomarkers,
                                   ontology = pcgr_data[["phenotype_ontology"]][["oncotree"]],
                                   primary_site = tumor_type)
    rlogging::message(paste0("Looking up SNV/InDel biomarkers for precision oncology - ", tumor_type))
  }

  ## get clinical evidence items that associated with query variants (non-regional - exact)
  var_eitems <- list()
  for(m in c('codon','exon','gene','exact')){
    var_eitems[[m]] <- data.frame()
  }
  for (db in c("civic","cgi")) {
    var_eitems_exact <- pcgrr::get_eitems_per_var(sample_calls, db = db,
                                           cancer_biomarkers = biomarkers,
                                            region_marker = F)
    var_eitems[['exact']] <- dplyr::bind_rows(var_eitems[['exact']], var_eitems_exact) %>%
      pcgrr::remove_cols_from_df(cnames = c("EITEM_CONSEQUENCE",
                                            "EITEM_CODON", "EITEM_EXON"))

  }

  ## get clinical evidence items that associated with query variants (regional - codon/exon/gene)
  eitems_region <- pcgrr::get_eitems_per_var(sample_calls, db = "civic",
                                          cancer_biomarkers = biomarkers,
                                          region_marker = T)

  ## for regional biomarkers - perform additional checks/filtering (making sure variants are of correct consequence, at the correct amino acid position etc)
  for(m in c('codon','exon','gene')){
    if(nrow(eitems_region) > 0){
      var_eitems[[m]] <- pcgrr::filter_eitems_per_var(eitems_region, marker_type = m)
    }
  }

  ## ignore variant biomarkers at the codon/exon/gene level if they are
  ## already present at the exact (variant level)
  if(nrow(var_eitems[['exact']]) > 0){
    for(m in c('codon','exon','gene')){
      if(nrow(var_eitems[[m]] > 0)){
        var_eitems[[m]] <- var_eitems[[m]] %>%
          dplyr::anti_join(dplyr::select(var_eitems[['exact']], GENOMIC_CHANGE),
                           by = c("GENOMIC_CHANGE"))
      }
    }
  }

  ## ignore variant biomarkers at the exon/gene level if they are
  ## already present at the codon level
  if(nrow(var_eitems[['codon']]) > 0){
    for(m in c('exon','gene')){
      if(nrow(var_eitems[[m]] > 0)){
        var_eitems[[m]] <- var_eitems[[m]] %>%
          dplyr::anti_join(dplyr::select(var_eitems[['codon']], GENOMIC_CHANGE),
                           by = c("GENOMIC_CHANGE"))
      }
    }
  }

  ## limit evidence items to exact/codon and exon (ignore biomarkers reported with a gene-level resolution)
  all_var_evidence_items <- all_var_evidence_items %>%
    dplyr::bind_rows(var_eitems[['exact']]) %>%
    dplyr::bind_rows(var_eitems[['codon']]) %>%
    dplyr::bind_rows(var_eitems[['exon']])


  ## log the types and number of clinical evidence items found (exact / codon / exon)
  rlogging::message(paste0("Found n = ",nrow(var_eitems[['exact']]), " clinical evidence item(s) at the variant level, ",
                           length(unique(var_eitems[['exact']]$GENOMIC_CHANGE))," unique variant(s)"))
  if(nrow(var_eitems[['exact']]) > 0){
    rlogging::message(paste0("Variants: ", paste(unique(paste(var_eitems[['exact']]$SYMBOL,
                                                            var_eitems[['exact']]$CONSEQUENCE,
                                                            var_eitems[['exact']]$PROTEIN_CHANGE,sep=":")),
                                                collapse=", ")))
  }

  rlogging::message(paste0("Found n = ",nrow(var_eitems[['codon']]), " other clinical evidence item(s) at the codon level, ",
                           length(unique(var_eitems[['codon']]$GENOMIC_CHANGE))," unique variant(s)"))
  if(nrow(var_eitems[['codon']]) > 0){
    rlogging::message(paste0("Variants: ", paste(unique(paste(var_eitems[['codon']]$SYMBOL,
                                                            var_eitems[['codon']]$CONSEQUENCE,
                                                            var_eitems[['codon']]$PROTEIN_CHANGE,sep=":")),
                                                collapse=", ")))
  }

  rlogging::message(paste0("Found n = ",nrow(var_eitems[['exon']]), " other clinical evidence item(s) at the exon level, ",
                           length(unique(var_eitems[['exon']]$GENOMIC_CHANGE))," unique variant(s)"))

  if(nrow(var_eitems[['exon']]) > 0){
    rlogging::message(paste0("Variants: ", paste(unique(paste(var_eitems[['exon']]$SYMBOL,
                                                              var_eitems[['exon']]$CONSEQUENCE,
                                                              var_eitems[['exon']]$PROTEIN_CHANGE,sep=":")),
                                                 collapse=", ")))
  }

  ## Organize all variants in a list object 'clin_items', organized through
  ## 1) tumor type (specific_ttype|any_ttype|other_ttype)
  ## 2) evidence type (diagnostic|prognostic|predictive)
  ## 3) clinical significance ('A_B','C_D_E','any')

  clin_eitems <-
    pcgrr::get_eitems_per_var_by_significance(all_var_evidence_items, pcgr_data = pcgr_data,
                                              alteration_type = "MUT")

  variant_set <- data.frame()
  if(nrow(all_var_evidence_items) > 0){
    variant_tags <-
      pcgr_data[["annotation_tags"]][["all"]][!pcgr_data[["annotation_tags"]][["all"]]
                                              %in% c("CIVIC_ID","CIVIC_ID_SEGMENT",
                                                     "CGI_ID_SEGMENT","CGI_ID")]
    variant_set <-  all_var_evidence_items %>%
      dplyr::select(dplyr::one_of(variant_tags)) %>%
      dplyr::distinct()
  }

  return(list("clin_eitem" = clin_eitems, "variant_set" = variant_set))

}


#' Function that retrieves clinical evidence items (CIVIC, CBMDB) for CNA aberrations
#'
#' @param onco_ts_sets data frame with lost tumor suppressor genes and gained oncogenes
#' @param pcgr_data object with PCGR annotation data
#' @param tumor_type Primary site of tumor
#' @param tumor_type_specificity run against specific tumor type ('specific'), as specified in pcgr_config or any tumor type 'any'
#'
#' @return list
#'
get_clin_assocs_cna <- function(onco_ts_sets, pcgr_data, tumor_type,
                                tumor_type_specificity = "any") {

  assertthat::assert_that("oncogene_gain" %in% names(onco_ts_sets) & "tsgene_loss" %in% names(onco_ts_sets),
                          msg = "Input object 'onco_ts_sets' does not contain appropriate lists ('oncogene_gain' and 'tsgene_loss')")
  invisible(assertthat::assert_that(tumor_type_specificity == "any" | tumor_type_specificity == "specific",
                                    msg = paste0("Argument 'tumor_type_specificy' can only take two values: 'any' or 'specific' and NOT: ",
                                                 tumor_type_specificity)))
  if(tumor_type_specificity == "specific"){
    invisible(assertthat::assert_that(tumor_type != "Cancer, NOS",
                                      msg = paste0("If tumor_specificy is 'specific', 'tumor_type' must be specified and NOT:  ",
                                                   tumor_type)))
  }

  variant_set <- data.frame()

  ## get all CNA clinical evidence items (civic, cgi)
  biomarkers <-
    pcgrr::load_all_eitems(pcgr_data = pcgr_data, alteration_type = "CNA", origin = "Somatic")

  ## limit CNA clinical evidence items by primary tumor site
  if (tumor_type_specificity == "any") {
    rlogging::message(paste0("Looking up sCNA biomarkers for precision oncology - any tumortype"))
  }else{
    biomarkers <-
      pcgrr::filter_eitems_by_site(biomarkers,
                                   ontology = pcgr_data[["phenotype_ontology"]][["oncotree"]],
                                   primary_site = tumor_type)
    rlogging::message(paste0("Looking up sCNA biomarkers for precision oncology - ", tumor_type))
  }

  ## intersect known CNA clinical evidence items with those found in the queryset
  for (type in c("tsgene_loss", "oncogene_gain")) {
    if ( NROW(onco_ts_sets[[type]]) > 0) {
      biomarker_hits <- as.data.frame(
        dplyr::inner_join(biomarkers,
                         onco_ts_sets[[type]],
                         by = c("SYMBOL", "CNA_TYPE"))
      )
      if (nrow(biomarker_hits) > 0) {
        variant_set <- dplyr::bind_rows(variant_set, biomarker_hits) %>%
          dplyr::select(dplyr::one_of(pcgr_data[['annotation_tags']][['cna_display']])) %>%
          dplyr::distinct()
      }
    }
  }


  ## Organize all variants in a list object 'clin_items', organized through
  ## 1) tumor type (specific_ttype|any_ttype|other_ttype)
  ## 2) evidence type (diagnostic|prognostic|predictive)
  ## 3) clinical significance ('A_B','C_D_E','any')

  clin_eitems <- pcgrr::get_eitems_per_var_by_significance(variant_set, pcgr_data = pcgr_data,
                                                             alteration_type = "CNA")


  return(list("clin_eitem" = clin_eitems,
              "variant_set" = variant_set))

}


load_all_eitems <- function(pcgr_data = NULL, alteration_type = "MUT", origin = "Somatic"){

  biomarkers <- list()
  for (db in c("civic", "cgi")) {

    invisible(assertthat::assert_that(!is.null(pcgr_data[['biomarkers']]), msg = "pcgr_data[['biomarkers']] is NULL - existing"))
    invisible(assertthat::assert_that(db %in% names(pcgr_data[['biomarkers']]),
                                      msg = paste0("Datasource ", db, " cannot be found in pcgr_data[['biomarkers']]")))
    invisible(assertthat::assert_that(is.data.frame(pcgr_data[['biomarkers']][[db]]),
                                      msg = paste0("Object pcgr_data[['biomarkers']][['", db, "']] must be of type data frame, not ",
                                      class(pcgr_data[['biomarkers']][[db]]))))

    if(alteration_type == 'CNA'){
      biomarkers[[db]] <-
        pcgr_data[["biomarkers"]][[db]] %>%
          dplyr::filter(ALTERATION_TYPE == alteration_type & !is.na(EITEM_CONSEQUENCE) &
                          stringr::str_detect(VARIANT_ORIGIN, origin)) %>%
          dplyr::rename(CNA_TYPE = EITEM_CONSEQUENCE) %>%
          pcgrr::remove_cols_from_df(cnames = c("VARIANT_NAME", "STATUS", "DRUG_INTERACTION_TYPE",
                                        "ALTERATION_TYPE", "EITEM_CODON", "EITEM_EXON",
                                        "MOLECULE_CHEMBL_ID", "MAPPING_RANK",
                                        "BIOMARKER_DESCRIPTION"))
    }
    if(alteration_type == 'MUT'){
      biomarkers[[db]] <-
        pcgr_data[["biomarkers"]][[db]] %>%
        dplyr::filter(ALTERATION_TYPE == alteration_type &
                        stringr::str_detect(VARIANT_ORIGIN, origin)) %>%
        pcgrr::remove_cols_from_df(cnames = c("VARIANT_NAME", "STATUS", "DRUG_INTERACTION_TYPE",
                                     "ALTERATION_TYPE", "MOLECULE_CHEMBL_ID",
                                     "MAPPING_RANK","BIOMARKER_DESCRIPTION"))
    }
    biomarkers[[db]] <- biomarkers[[db]] %>%
      dplyr::mutate(SOURCE_DB = db)

  }

  all_biomarkers <- dplyr::bind_rows(biomarkers[['civic']], biomarkers[['cgi']])


  return(all_biomarkers)
}


get_eitems_per_var <- function(sample_calls, db = "civic",
                               cancer_biomarkers = NULL, region_marker = T){

  invisible(assertthat::assert_that(db == "civic" | db == "cgi", msg = "Argument 'db' can be one of civic or cgi"))
  invisible(assertthat::assert_that(!is.null(cancer_biomarkers), msg = "Argument 'cancer_biomarkers' cannot be NULL"))
  invisible(assertthat::assert_that(is.data.frame(cancer_biomarkers),
                                    msg = paste0("Argument 'cancer_biomarkers' must be of type data frame, not ",
                                    class(cancer_biomarkers))))
  assertable::assert_colnames(cancer_biomarkers, c("EVIDENCE_ID","SYMBOL"), only_colnames = F, quiet = T)

  eitems <- data.frame()

  if(region_marker == F){
    if(db == 'civic'){
      assertable::assert_colnames(sample_calls, c("CIVIC_ID","CIVIC_ID_SEGMENT","SYMBOL"),
                                  only_colnames = F, quiet = T)
      sample_calls_db <- sample_calls %>% dplyr::filter(!is.na(CIVIC_ID))
      if (nrow(sample_calls_db) > 0) {
        eitems <- sample_calls_db %>%
          tidyr::separate_rows(CIVIC_ID, sep = ",") %>%
          dplyr::select(dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
          dplyr::rename(EVIDENCE_ID = CIVIC_ID) %>%
          dplyr::mutate(EVIDENCE_ID = as.character(EVIDENCE_ID)) %>%
          pcgrr::remove_cols_from_df(cnames = c('CIVIC_ID_SEGMENT'))
      }
    }else{
      assertable::assert_colnames(sample_calls, c("CGI_ID","CGI_ID_SEGMENT","SYMBOL"),
                                  only_colnames = F, quiet = T)
      sample_calls_db <- sample_calls %>% dplyr::filter(!is.na(CGI_ID))
      if (nrow(sample_calls_db) > 0) {
        eitems <- sample_calls_db %>%
          tidyr::separate_rows(CGI_ID, sep = ",") %>%
          dplyr::select(dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
          dplyr::rename(EVIDENCE_ID = CGI_ID) %>%
          dplyr::mutate(EVIDENCE_ID = as.character(EVIDENCE_ID)) %>%
          pcgrr::remove_cols_from_df(cnames = c('CGI_ID_SEGMENT'))
      }
    }
  }else{
    if(db == 'civic'){
      assertable::assert_colnames(sample_calls, c("CIVIC_ID","CIVIC_ID_SEGMENT","SYMBOL"),
                                  only_colnames = F, quiet = T)
      sample_calls_db <- sample_calls %>% dplyr::filter(!is.na(CIVIC_ID_SEGMENT))
      if (nrow(sample_calls_db) > 0) {
        eitems <- sample_calls_db %>%
          tidyr::separate_rows(CIVIC_ID_SEGMENT, sep = ",") %>%
          dplyr::select(dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
          dplyr::rename(EVIDENCE_ID = CIVIC_ID_SEGMENT) %>%
          dplyr::mutate(EVIDENCE_ID = as.character(EVIDENCE_ID)) %>%
          pcgrr::remove_cols_from_df(cnames = c('CIVIC_ID'))
      }
    }else{
      assertable::assert_colnames(sample_calls, c("CGI_ID","CGI_ID_SEGMENT","SYMBOL"),
                                  only_colnames = F, quiet = T)
      sample_calls_db <- sample_calls %>% dplyr::filter(!is.na(CGI_ID_SEGMENT))
      if (nrow(sample_calls_db) > 0) {
        eitems <- sample_calls_db %>%
          tidyr::separate_rows(CGI_ID_SEGMENT, sep = ",") %>%
          dplyr::select(dplyr::one_of(pcgr_data[["annotation_tags"]][["all"]])) %>%
          dplyr::rename(EVIDENCE_ID = CGI_ID_SEGMENT) %>%
          dplyr::mutate(EVIDENCE_ID = as.character(EVIDENCE_ID)) %>%
          pcgrr::remove_cols_from_df(cnames = c('CGI_ID'))
      }
    }
  }

  if(nrow(eitems) > 0){
    eitems <- eitems %>%
      dplyr::inner_join(cancer_biomarkers, by = c("EVIDENCE_ID", "SYMBOL")) %>%
      dplyr::distinct()

    if(db == 'civic'){
      eitems <- eitems %>%
        pcgrr::remove_cols_from_df(cnames = c("CGI_ID","CGI_ID_SEGMENT"))
    }

    if(db == 'cgi'){
      eitems <- eitems %>%
        pcgrr::remove_cols_from_df(cnames = c("CIVIC_ID","CIVIC_ID_SEGMENT"))
    }
  }


  return(eitems)

}

filter_eitems_per_var <- function(eitems, marker_type = "codon"){

  invisible(assertthat::assert_that(marker_type == "codon" | marker_type == "exon" | marker_type == "gene",
                          msg = "Argument marker_type can only be any of 'exon','codon' or 'gene'"))
  invisible(assertthat::assert_that(is.data.frame(eitems), msg = "Argument eitems must be of type data.frame()"))
  assertable::assert_colnames(eitems, c("EITEM_CODON", "EITEM_CONSEQUENCE", "AMINO_ACID_START",
                                        "BIOMARKER_MAPPING", "AMINO_ACID_END", "CONSEQUENCE",
                                        "CODING_STATUS", "EITEM_EXON", "SYMBOL"),
                              only_colnames = F,
                              quiet = T)

  filtered_eitems <- data.frame()
  if (marker_type == "codon") {
    if (nrow(eitems[!is.na(eitems$EITEM_CODON) & eitems$BIOMARKER_MAPPING == "codon", ]) > 0) {
      filtered_eitems <- eitems %>%
        dplyr::filter(!is.na(EITEM_CODON) & BIOMARKER_MAPPING == "codon") %>%
        dplyr::filter(EITEM_CODON <= AMINO_ACID_END & EITEM_CODON >= AMINO_ACID_START &
                        (!is.na(EITEM_CONSEQUENCE) & startsWith(CONSEQUENCE, EITEM_CONSEQUENCE) |
                           is.na(EITEM_CONSEQUENCE)) & CODING_STATUS == "coding")
    }
  }

  if (marker_type == "exon") {
    if (nrow(eitems[!is.na(eitems$EITEM_EXON) & eitems$BIOMARKER_MAPPING == "exon", ]) > 0) {
      filtered_eitems <- eitems %>%
        dplyr::filter(!is.na(EITEM_EXON) & BIOMARKER_MAPPING == "exon") %>%
        dplyr::filter(EXON == EITEM_EXON &
                        (!is.na(EITEM_CONSEQUENCE) &
                          startsWith(CONSEQUENCE, EITEM_CONSEQUENCE) |
                         is.na(EITEM_CONSEQUENCE)) & CODING_STATUS == "coding")
    }
  }

  if (marker_type == "gene") {
    if (nrow(eitems[eitems$BIOMARKER_MAPPING == "gene", ]) > 0) {
      filtered_eitems <- eitems %>%
        dplyr::filter(BIOMARKER_MAPPING == "gene" & CODING_STATUS == "coding")
    }
  }

  if (nrow(filtered_eitems) > 0) {
    filtered_eitems <- filtered_eitems %>%
      pcgrr::remove_cols_from_df(cnames = c("EITEM_CONSEQUENCE",
                                            "EITEM_CODON", "EITEM_EXON"))
      #dplyr::select(-c(EITEM_EXON, EITEM_CODON, EITEM_CONSEQUENCE))
  }
  return(filtered_eitems)

}

filter_eitems_by_site <- function(biomarkers = NULL, ontology = NULL, primary_site = "") {

  invisible(assertthat::assert_that(nchar(primary_site) > 0, msg = "Argument 'primary_site' cannot be an empty string"))
  invisible(assertthat::assert_that(!is.null(biomarkers), msg = "Argument 'biomarkers' cannot be NULL"))
  invisible(assertthat::assert_that(is.data.frame(biomarkers),
                                    msg = paste0("Argument 'biomarkers' must be of type data frame, not ",
                                                 class(biomarkers))))
  invisible(assertthat::assert_that(!is.null(ontology), msg = "Argument 'biomarkers' cannot be NULL"))
  invisible(assertthat::assert_that(is.data.frame(ontology),
                                    msg = paste0("Argument 'ontology' must be of type data frame, not ",
                                                 class(ontology))))
  assertable::assert_colnames(biomarkers, c("DISEASE_ONTOLOGY_ID"), only_colnames = F, quiet = T)
  assertable::assert_colnames(ontology, c("primary_site","do_id","cui","cui_name"), only_colnames = F, quiet = T)

  tumor_phenotypes_site <-
    dplyr::semi_join(dplyr::select(ontology,
                                   primary_site, do_id, cui, cui_name),
                     data.frame("primary_site" = primary_site, stringsAsFactors = F),
                     by = c("primary_site" = "primary_site")) %>%
    dplyr::filter(!is.na(do_id)) %>%
    dplyr::distinct()

  biomarkers <- biomarkers %>%
    dplyr::semi_join(tumor_phenotypes_site, by = c("DISEASE_ONTOLOGY_ID" = "do_id"))

  return(biomarkers)
}

get_eitems_per_var_by_significance <- function(var_eitems, pcgr_data, alteration_type = "MUT"){

  clin_eitems_list <- list()
  for (type in c("diagnostic", "prognostic", "predictive")) {
    clin_eitems_list[[type]] <- list()
    for (elevel in c("any", "A_B", "C_D_E")) {
      clin_eitems_list[[type]][[elevel]] <- data.frame()
    }
  }

  tags_display <- pcgr_data[["annotation_tags"]][["tier_1_2_display"]]
  if(alteration_type == 'CNA'){
    tags_display <- pcgr_data[['annotation_tags']][['cna_display']]
  }

  if (nrow(var_eitems) > 0) {
    for (type in c("prognostic", "diagnostic", "predictive")) {
      clin_eitems_list[[type]][["any"]] <- var_eitems %>%
        dplyr::select(dplyr::one_of(tags_display)) %>%
        dplyr::filter(EVIDENCE_TYPE == stringr::str_to_title(type)) %>%
        dplyr::arrange(EVIDENCE_LEVEL, desc(RATING))
      if (nrow(clin_eitems_list[[type]][["any"]]) > 0) {
        if (nrow(clin_eitems_list[[type]][["any"]][stringr::str_detect(clin_eitems_list[[type]][["any"]]$EVIDENCE_LEVEL, "^(A|B|B1|B2):"),]) > 0) {
          clin_eitems_list[[type]][["A_B"]] <- clin_eitems_list[[type]][["any"]] %>%
            dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL, "^(A|B|B1|B2):")) %>%
            dplyr::arrange(EVIDENCE_LEVEL, desc(RATING))
        }
        if (nrow(clin_eitems_list[[type]][["any"]][stringr::str_detect(clin_eitems_list[[type]][["any"]]$EVIDENCE_LEVEL, "^(C|D|E):"),]) > 0) {
          clin_eitems_list[[type]][["C_D_E"]] <- clin_eitems_list[[type]][["any"]] %>%
            dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL, "^(C|D|E):")) %>%
            dplyr::arrange(EVIDENCE_LEVEL, desc(RATING))
        }
      }
    }
  }
  return(clin_eitems_list)


}
