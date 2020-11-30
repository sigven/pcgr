#' Function that matches clinical evidence items (CIVIC, CBMDB)
#' against somatic cancer variants detected in tumor
#'
#' @param sample_calls data frame with sample variants
#' @param annotation_tags list with annotation tags for display in report
#' @param eitems data frame with clinical evidence items
#'
#' @return list

get_clin_assocs_snv_indel <- function(sample_calls,
                                      annotation_tags = NULL,
                                      eitems = NULL) {
                                      #tumor_type = NA,
                                      #tumor_type_specificity = "any",
                                      #mapping_stringency = 1) {

  invisible(assertthat::assert_that(!is.null(annotation_tags)))
  invisible(assertthat::assert_that(!is.null(eitems)))
  invisible(assertthat::assert_that(is.data.frame(sample_calls)))
  invisible(assertthat::assert_that(is.data.frame(eitems)))
  invisible(assertthat::assert_that(
    "all" %in% names(annotation_tags) &
      "tier_1_2_display" %in% names(annotation_tags)))

  ## Initialize lists that hold clinical evidence items -
  ## across evidence types (diagnostic/prognostic/predictive)
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

  var_eitems <- list()
  for (m in c("codon", "exon", "gene", "exact")) {
    var_eitems[[m]] <- data.frame()
  }

  ## get clinical evidence items that associated with
  ## query variants (non-regional - exact), civic + cgi
  for (db in c("civic", "cgi")) {
    var_eitems_exact <-
      pcgrr::match_eitems_to_var(
        sample_calls,
        db = db,
        colset = annotation_tags$all,
        eitems = eitems,
        region_marker = F)

    var_eitems[["exact"]] <- var_eitems[["exact"]] %>%
      dplyr::bind_rows(var_eitems_exact) %>%
      pcgrr::remove_cols_from_df(
        cnames = c("EITEM_CONSEQUENCE", "EITEM_CODON",
                   "EITEM_EXON"))

  }

  ## get clinical evidence items that associated with
  ## query variants (regional - codon/exon/gene)
  var_eitems_regional <-
    pcgrr::match_eitems_to_var(
      sample_calls,
      db = "civic",
      colset = annotation_tags$all,
      eitems = eitems,
      region_marker = T)

  ## for regional biomarkers - perform additional quality checks
  ## (making sure variants are of correct consequence,
  ## at the correct amino acid position etc)
  for (m in c("codon", "exon", "gene")) {
    if (nrow(var_eitems_regional) > 0) {
      var_eitems[[m]] <-
        pcgrr::qc_var_eitems(var_eitems = var_eitems_regional,
                             marker_type = m)
    }
  }


  var_eitems <- pcgrr::deduplicate_eitems(var_eitems = var_eitems,
                                          target_type = "exact",
                                          target_other =
                                            c("codon", "exon", "gene"))

  var_eitems <- pcgrr::deduplicate_eitems(var_eitems = var_eitems,
                                          target_type = "codon",
                                          target_other =
                                            c("exon", "gene"))

  ## limit evidence items to exact/codon and exon
  ## (ignore biomarkers reported with a gene-level resolution)
  all_var_evidence_items <- all_var_evidence_items %>%
    dplyr::bind_rows(var_eitems[["exact"]]) %>%
    dplyr::bind_rows(var_eitems[["codon"]]) %>%
    dplyr::bind_rows(var_eitems[["exon"]])

  ## log the types and number of clinical
  ## evidence items found (exact / codon / exon)
  pcgrr::log_var_eitem_stats(var_eitems = var_eitems, target_type = "exact")
  pcgrr::log_var_eitem_stats(var_eitems = var_eitems, target_type = "codon")
  pcgrr::log_var_eitem_stats(var_eitems = var_eitems, target_type = "exon")

  ## Organize all variants in a list object 'clin_items', organized through
  ## 1) tumor type (specific_ttype|any_ttype|other_ttype)
  ## 2) evidence type (diagnostic|prognostic|predictive)
  ## 3) clinical significance ('A_B','C_D_E','any')

  clin_eitems <-
    pcgrr::structure_var_eitems(
      all_var_evidence_items,
      annotation_tags = annotation_tags,
      alteration_type = "MUT")

  variant_set <- data.frame()
  if (nrow(all_var_evidence_items) > 0) {
    variant_tags <-
      annotation_tags[["all"]][!annotation_tags[["all"]]
                                              %in% c("CIVIC_ID",
                                                     "CIVIC_ID_SEGMENT",
                                                     "CGI_ID_SEGMENT",
                                                     "CGI_ID")]
    variant_set <-  all_var_evidence_items %>%
      dplyr::select(dplyr::one_of(variant_tags)) %>%
      dplyr::distinct()
  }

  return(list("clin_eitem" = clin_eitems, "variant_set" = variant_set))

}


#' Function that retrieves clinical evidence items (CIVIC, CBMDB) for
#' CNA aberrations
#'
#' @param onco_ts_sets data frame with annotations with respect to
#' lost tumor suppressor genes and gained oncogenes
#' @param annotation_tags list with annotation tags for display in report
#' @param eitems Data frame with clinical evidence items
#'
#' @return list
#'
get_clin_assocs_cna <- function(onco_ts_sets,
                                annotation_tags = NULL,
                                eitems = NULL){

  assertthat::assert_that(
    "oncogene_gain" %in% names(onco_ts_sets) &
      "tsgene_loss" %in% names(onco_ts_sets),
    msg = paste0(
      "Input object 'onco_ts_sets' does not contain ",
      "appropriate lists ('oncogene_gain' and 'tsgene_loss')"))

  invisible(assertthat::assert_that(!is.null(annotation_tags)))
  invisible(assertthat::assert_that(!is.null(eitems)))
  invisible(assertthat::assert_that(is.data.frame(eitems)))
  invisible(assertthat::assert_that(
    "all" %in% names(annotation_tags) &
      "cna_display" %in% names(annotation_tags)))

  assertable::assert_colnames(eitems,
                              colnames = c("SYMBOL", "CNA_TYPE"),
                              quiet = T, only_colnames = F)

  variant_set <- data.frame()

  ## intersect known CNA clinical evidence items with
  ## those found in the queryset
  for (type in c("tsgene_loss", "oncogene_gain")) {
    if (NROW(onco_ts_sets[[type]]) > 0) {

      assertable::assert_colnames(onco_ts_sets[[type]],
                                  colnames = c("SYMBOL", "CNA_TYPE"),
                                  quiet = T, only_colnames = F)

      eitem_hits <- as.data.frame(
        dplyr::inner_join(eitems,
                         onco_ts_sets[[type]],
                         by = c("SYMBOL", "CNA_TYPE"))
      )
      if (nrow(eitem_hits) > 0) {
        variant_set <- dplyr::bind_rows(variant_set, eitem_hits) %>%
          dplyr::select(
            dplyr::one_of(annotation_tags[["cna_display"]])) %>%
          dplyr::distinct()
      }
    }
  }


  ## Organize all variants in a list object 'clin_items', organized through
  ## 1) tumor type (specific_ttype|any_ttype|other_ttype)
  ## 2) evidence type (diagnostic|prognostic|predictive)
  ## 3) clinical significance ('A_B','C_D_E','any')

  clin_eitems <- pcgrr::structure_var_eitems(
    variant_set,
    annotation_tags = annotation_tags,
    alteration_type = "CNA")

  return(list("clin_eitem" = clin_eitems,
              "variant_set" = variant_set))

}

load_eitems <- function(eitems_raw = NULL,
                        ontology = NULL,
                        alteration_type = "MUT",
                        origin = "Somatic",
                        tumor_type_specificity = NULL,
                        tumor_type = NULL) {

  invisible(assertthat::assert_that(
    !is.null(eitems_raw),
    msg = "eitems_raw is NULL - existing"))

  invisible(
    assertthat::assert_that(
      origin == "Somatic" | origin == "Germline",
      msg = paste0("Argument 'origin' can only take ",
                   "two values: 'Germline' or 'Somatic' and NOT: ",
                   origin)))

  if(origin == "Somatic"){
    invisible(
      assertthat::assert_that(
        !is.null(tumor_type_specificity),
        msg = paste0(
          "If origin is 'Somatic', ",
          "'tumor_type_specificity' must not be NULL  ",
          tumor_type)))
  }

  invisible(
    assertthat::assert_that(
      alteration_type == "MUT" |
        alteration_type == "CNA" |
        alteration_type == "MUT_LOF",
      msg = paste0("Argument 'alteration_type' can only take ",
                   "two values: 'MUT' or 'CNA' or 'MUT_LOF', and NOT: ",
                   alteration_type)))

  invisible(
    assertthat::assert_that(
      tumor_type_specificity == "any" |
        tumor_type_specificity == "specific",
      msg = paste0("Argument 'tumor_type_specificy' can only take ",
                   "two values: 'any' or 'specific' and NOT: ",
                   tumor_type_specificity)))
  if (tumor_type_specificity == "specific") {
    invisible(
      assertthat::assert_that(
        !is.na(tumor_type) & tumor_type != "Cancer, NOS",
        msg = paste0(
          "If tumor_type_specificity is 'specific', ",
          "'tumor_type' must be specified and NOT:  ",
          tumor_type)))
  }

  ## load all clinical evidence items (civic and cgi), by
  ## mutation type and origin
  eitems <-
    pcgrr::load_all_eitems(
      eitems = eitems_raw,
      alteration_type = alteration_type,
      origin = origin)


  if (tumor_type_specificity == "any") {
    rlogging::message(
      paste0(
        "Loading ", alteration_type, " biomarkers for precision oncology",
        "- any tumortype"))
  }else{

    ## limit clinical evidence items by primary tumor site if this
    ## is specified in arguments
    invisible(
      assertthat::assert_that(
        !is.null(ontology),
        msg = "Argument 'ontology' cannot be NULL"))
    invisible(assertthat::assert_that(
      is.data.frame(ontology),
      msg = paste0("Argument 'ontology' must be of type data frame, not ",
                   class(ontology))))
    assertable::assert_colnames(
      eitems, c("DISEASE_ONTOLOGY_ID"),
      only_colnames = F, quiet = T)
    assertable::assert_colnames(
      ontology,
      c("primary_site", "do_id", "cui", "cui_name"),
      only_colnames = F, quiet = T)

    eitems <-
      pcgrr::filter_eitems_by_site(
        eitems,
        ontology = ontology,
        primary_site = tumor_type)
    rlogging::message(
      paste0("Loading ", alteration_type,
             " for precision oncology - ",
             tumor_type))
  }
  return(eitems)


}

load_all_eitems <- function(eitems_raw = NULL,
                            alteration_type = "MUT",
                            origin = "Somatic") {

  invisible(
    assertthat::assert_that(
      origin == "Somatic" | origin == "Germline",
      msg = paste0("Argument 'origin' can only take ",
                   "two values: 'Germline' or 'Somatic' and NOT: ",
                   origin)))
  invisible(
    assertthat::assert_that(
      alteration_type == "MUT" |
        alteration_type == "CNA" |
        alteration_type == "MUT_LOF",
      msg = paste0("Argument 'alteration_type' can only take ",
                   "two values: 'MUT' or 'CNA' or 'MUT_LOF' and NOT: ",
                   alteration_type)))

  invisible(assertthat::assert_that(
    !is.null(eitems_raw),
    msg = "eitems_raw is NULL - existing"))

  selected_eitems <- list()
  for (db in c("civic", "cgi")) {

    invisible(assertthat::assert_that(
      db %in% names(eitems_raw),
      msg = paste0("Datasource ", db,
                   " cannot be found in eitems_raw (i.e. eitems_raw$db is NULL)")))
    invisible(assertthat::assert_that(
      is.data.frame(eitems_raw[[db]]),
      msg = paste0("Object eitems_raw[['",
                   db, "']] must be of type data frame, not ",
                   class(eitems_raw[[db]]))))

    assertable::assert_colnames(
      eitems_raw[[db]],
      c('ALTERATION_TYPE',
        'EITEM_CONSEQUENCE',
        'VARIANT_ORIGIN'),
      only_colnames = F,
      quiet = T)

    if(alteration_type == "CNA") {
      selected_eitems[[db]] <-
        eitems_raw[[db]] %>%
          dplyr::filter(ALTERATION_TYPE == alteration_type &
                          !is.na(EITEM_CONSEQUENCE) &
                          stringr::str_detect(VARIANT_ORIGIN, origin)) %>%
          dplyr::rename(CNA_TYPE = EITEM_CONSEQUENCE) %>%
          pcgrr::remove_cols_from_df(
            cnames =
              c("VARIANT_NAME",
                "STATUS",
                "DRUG_INTERACTION_TYPE",
                "EITEM_CODON",
                "EITEM_EXON",
                "MOLECULE_CHEMBL_ID",
                "MAPPING_RANK",
                "BIOMARKER_DESCRIPTION"))
    }
    if (alteration_type == "MUT" | alteration_type == "MUT_LOF") {
      selected_eitems[[db]] <-
        eitems_raw[[db]] %>%
        dplyr::filter(ALTERATION_TYPE == alteration_type &
                        stringr::str_detect(VARIANT_ORIGIN, origin)) %>%
        pcgrr::remove_cols_from_df(
          cnames =
            c("VARIANT_NAME",
              "STATUS",
              "DRUG_INTERACTION_TYPE",
              "MOLECULE_CHEMBL_ID",
              "MAPPING_RANK",
              "BIOMARKER_DESCRIPTION"))
    }
    selected_eitems[[db]] <- selected_eitems[[db]] %>%
      dplyr::mutate(SOURCE_DB = db)

  }

  all_eitems <- dplyr::bind_rows(selected_eitems[["civic"]],
                                 selected_eitems[["cgi"]])

  return(all_eitems)
}


match_eitems_to_var <- function(sample_calls,
                               db = "civic",
                               colset = NULL,
                               eitems = NULL,
                               region_marker = T) {

  invisible(assertthat::assert_that(
    db == "civic" | db == "cgi",
    msg = "Argument 'db' can be one of 'civic' or 'cgi'"))
  invisible(assertthat::assert_that(
    !is.null(eitems),
    msg = "Argument 'eitems' cannot be NULL"))
  invisible(assertthat::assert_that(
    is.data.frame(eitems),
    msg = paste0(
      "Argument 'eitems' must be of type data frame, not ",
      class(eitems))))
  invisible(assertthat::assert_that(
    !is.null(sample_calls) & is.data.frame(sample_calls),
    msg = paste0(
      "Argument 'sample_calls' must be of type data frame, not ",
      class(sample_calls))))
  assertable::assert_colnames(
    eitems, c("EVIDENCE_ID", "SYMBOL"),
    only_colnames = F, quiet = T)

  invisible(assertthat::assert_that(!is.null(colset)))
  invisible(assertthat::assert_that(is.character(colset)))

  # invisible(assertthat::assert_that("all" %in%
  #                                     names(annotation_tags)))

  evidence_identifiers <- c("CIVIC_ID", "CIVIC_ID_SEGMENT")
  if (region_marker == T) {
    evidence_identifiers <- c("CIVIC_ID_SEGMENT", "CIVIC_ID")
  }
  if(db == "cgi"){
    evidence_identifiers <- c("CGI_ID", "CGI_ID_SEGMENT")
    if (region_marker == T) {
      evidence_identifiers <- c("CGI_ID_SEGMENT", "CGI_ID")
    }
  }

  var_eitems <- data.frame()

  assertable::assert_colnames(
    sample_calls,
    c(evidence_identifiers, "SYMBOL"),
    only_colnames = F, quiet = T)

  sample_calls_db <- sample_calls %>%
    dplyr::filter(!is.na(!!sym(evidence_identifiers[1])))
  if (nrow(sample_calls_db) > 0) {
    var_eitems <- as.data.frame(sample_calls_db %>%
      tidyr::separate_rows(!!sym(evidence_identifiers[1]), sep = ",") %>%
      dplyr::select(
        dplyr::one_of(colset)) %>%
      dplyr::rename(EVIDENCE_ID = !!sym(evidence_identifiers[1])) %>%
      dplyr::mutate(EVIDENCE_ID = as.character(EVIDENCE_ID)) %>%
      pcgrr::remove_cols_from_df(cnames = evidence_identifiers[2])
    )
  }


  if (nrow(var_eitems) > 0) {
    var_eitems <- as.data.frame(var_eitems %>%
      dplyr::inner_join(eitems,
                        by = c("EVIDENCE_ID", "SYMBOL")) %>%
      dplyr::distinct() %>%
      pcgrr::remove_cols_from_df(cnames = evidence_identifiers)
    )
  }


  return(var_eitems)

}

##
##
qc_var_eitems <- function(var_eitems = NULL,
                          marker_type = "codon") {

  invisible(assertthat::assert_that(!is.null(var_eitems)))
  invisible(
    assertthat::assert_that(
      marker_type == "codon" |
        marker_type == "exon" |
        marker_type == "gene",
      msg = "Argument marker_type can only be any of 'exon','codon' or 'gene'"))
  invisible(
    assertthat::assert_that(
      is.data.frame(var_eitems),
      msg = "Argument eitems must be of type data.frame()"))
  assertable::assert_colnames(
    var_eitems, c("EITEM_CODON", "EITEM_CONSEQUENCE", "AMINO_ACID_START",
              "BIOMARKER_MAPPING", "AMINO_ACID_END", "CONSEQUENCE",
              "CODING_STATUS", "EITEM_EXON", "SYMBOL","EXON"),
    only_colnames = F,
    quiet = T)

  filtered_var_eitems <- data.frame()
  if (marker_type == "codon") {
    if (nrow(var_eitems[!is.na(var_eitems$EITEM_CODON) &
                    var_eitems$BIOMARKER_MAPPING == "codon", ]) > 0) {
      filtered_var_eitems <- var_eitems %>%
        dplyr::filter(!is.na(EITEM_CODON) & BIOMARKER_MAPPING == "codon") %>%
        dplyr::filter(EITEM_CODON <= AMINO_ACID_END &
                        EITEM_CODON >= AMINO_ACID_START &
                        (!is.na(EITEM_CONSEQUENCE) &
                           startsWith(CONSEQUENCE, EITEM_CONSEQUENCE) |
                           is.na(EITEM_CONSEQUENCE)) &
                        CODING_STATUS == "coding")

    }
  }

  if (marker_type == "exon") {
    if (nrow(var_eitems[!is.na(var_eitems$EITEM_EXON) &
                    var_eitems$BIOMARKER_MAPPING == "exon", ]) > 0) {
      filtered_var_eitems <- var_eitems %>%
        dplyr::filter(!is.na(EITEM_EXON) & BIOMARKER_MAPPING == "exon") %>%
        dplyr::filter(EXON == EITEM_EXON &
                        (!is.na(EITEM_CONSEQUENCE) &
                          startsWith(CONSEQUENCE, EITEM_CONSEQUENCE) |
                         is.na(EITEM_CONSEQUENCE)) & CODING_STATUS == "coding")
    }
  }

  if (marker_type == "gene") {
    if (nrow(var_eitems[var_eitems$BIOMARKER_MAPPING == "gene", ]) > 0) {
      filtered_var_eitems <- var_eitems %>%
        dplyr::filter(BIOMARKER_MAPPING == "gene" &
                        CODING_STATUS == "coding") %>%
        dplyr::filter((!is.na(EITEM_CONSEQUENCE) &
                         startsWith(CONSEQUENCE, EITEM_CONSEQUENCE) |
                         is.na(EITEM_CONSEQUENCE)) & CODING_STATUS == "coding")

    }
  }

  if (nrow(filtered_var_eitems) > 0) {

    if("LOSS_OF_FUNCTION" %in% colnames(filtered_var_eitems) &
       "ALTERATION_TYPE" %in% colnames(filtered_var_eitems)){

      filtered_var_eitems <- filtered_var_eitems %>%
        dplyr::filter((LOSS_OF_FUNCTION == T &
                         ALTERATION_TYPE == "MUT_LOF") |
                        is.na(LOSS_OF_FUNCTION) |
                        (LOSS_OF_FUNCTION == F &
                           ALTERATION_TYPE != "MUT_LOF"))
    }
  }

  filtered_var_eitems <- filtered_var_eitems %>%
    pcgrr::remove_cols_from_df(cnames = c("EITEM_CONSEQUENCE",
                                          "EITEM_CODON",
                                          "EITEM_EXON"))

  return(filtered_var_eitems)

}

filter_eitems_by_site <- function(eitems = NULL,
                                  ontology = NULL,
                                  primary_site = "") {

  invisible(
    assertthat::assert_that(
      nchar(primary_site) > 0,
      msg = "Argument 'primary_site' cannot be an empty string"))
  invisible(
    assertthat::assert_that(
      !is.null(eitems),
      msg = "Argument 'eitems' cannot be NULL"))
  invisible(
    assertthat::assert_that(
      is.data.frame(eitems),
      msg = paste0("Argument 'eitems' must be of type data frame, not ",
                   class(eitems))))
  invisible(
    assertthat::assert_that(
      !is.null(ontology),
      msg = "Argument 'ontology' cannot be NULL"))
  invisible(assertthat::assert_that(
    is.data.frame(ontology),
    msg = paste0("Argument 'ontology' must be of type data frame, not ",
                 class(ontology))))
  assertable::assert_colnames(
    eitems, c("DISEASE_ONTOLOGY_ID"),
    only_colnames = F, quiet = T)
  assertable::assert_colnames(
    ontology,
    c("primary_site", "do_id", "cui", "cui_name"),
    only_colnames = F, quiet = T)
  invisible(assertthat::assert_that(
    primary_site %in% unique(ontology$primary_site),
    msg = paste0("Tumor primary site ", primary_site,
                 " is not a recognized primary site",
                 " (possible values: ", paste(unique(ontology$primary_site),
                                              collapse=", "),")")))

  tumor_phenotypes_site <-
    dplyr::semi_join(
      dplyr::select(ontology,
                    primary_site, do_id, cui, cui_name),
      data.frame("primary_site" = primary_site, stringsAsFactors = F),
      by = c("primary_site" = "primary_site")) %>%
    dplyr::filter(!is.na(do_id)) %>%
    dplyr::distinct()

  eitems <- eitems %>%
    dplyr::semi_join(tumor_phenotypes_site,
                     by = c("DISEASE_ONTOLOGY_ID" = "do_id"))

  return(eitems)
}

structure_var_eitems <- function(var_eitems,
                              annotation_tags,
                              alteration_type = "MUT") {

  clin_eitems_list <- list()
  for (type in c("diagnostic", "prognostic", "predictive")) {
    clin_eitems_list[[type]] <- list()
    for (elevel in c("any", "A_B", "C_D_E")) {
      clin_eitems_list[[type]][[elevel]] <- data.frame()
    }
  }

  tags_display <- annotation_tags[["tier_1_2_display"]]
  if (alteration_type == "CNA") {
    tags_display <- annotation_tags[["cna_display"]]
  }

  if (nrow(var_eitems) > 0) {
    for (type in c("prognostic", "diagnostic", "predictive")) {
      clin_eitems_list[[type]][["any"]] <- var_eitems %>%
        dplyr::select(dplyr::one_of(tags_display)) %>%
        dplyr::filter(EVIDENCE_TYPE == stringr::str_to_title(type)) %>%
        dplyr::arrange(EVIDENCE_LEVEL, desc(RATING))
      if (nrow(clin_eitems_list[[type]][["any"]]) > 0) {
        clin_eitems_list[[type]][["A_B"]] <-
          clin_eitems_list[[type]][["any"]] %>%
          dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL, "^(A|B|B1|B2):"))

        if (NROW(clin_eitems_list[[type]][["A_B"]]) > 0) {
          clin_eitems_list[[type]][["A_B"]] <-
            clin_eitems_list[[type]][["A_B"]] %>%
            dplyr::arrange(EVIDENCE_LEVEL, desc(RATING))
        }

        clin_eitems_list[[type]][["C_D_E"]] <-
          clin_eitems_list[[type]][["any"]] %>%
          dplyr::filter(stringr::str_detect(EVIDENCE_LEVEL, "^(C|D|E):"))

        if (NROW(clin_eitems_list[[type]][["C_D_E"]]) > 0) {
          clin_eitems_list[[type]][["C_D_E"]] <-
            clin_eitems_list[[type]][["C_D_E"]] %>%
            dplyr::arrange(EVIDENCE_LEVEL, desc(RATING))
        }
      }
    }
  }
  return(clin_eitems_list)


}

deduplicate_eitems <- function(var_eitems = NULL,
                               target_type = "exact",
                               target_other = c("codon","exon","gene")){

  invisible(
    assertthat::assert_that(!is.null(var_eitems),
                          msg = "Object 'var_eitems' cannot be NULL"))
  invisible(
    assertthat::assert_that(target_type == "exact" | target_type == "codon",
                            msg = paste0("Argument 'target_type' can only",
                                         "take on values 'codon' or 'exact'")))

  if(target_type == "exact"){
    invisible(
      assertthat::assert_that(
        ("codon" %in% target_other &
           "exon" %in% target_other &
           "gene" %in% target_other &
           length(target_other) == 3),
        msg = paste0("Argument target_other must be ",
                     "specified as c('codon','exon','gene')"))
    )
  }else{
    invisible(
      assertthat::assert_that(
        ("exon" %in% target_other &
           "gene" %in% target_other &
           length(target_other) == 2),
        msg = paste0("Argument target_other must be ",
                     "specified as c('exon','gene')"))
    )
  }


  ## ignore variant evidence items at the codon/exon/gene level if they are
  ## already present at the exact (variant level)
  ##
  ## OR
  ##
  ## ignore variant biomarkers at the exon/gene level if they are
  ## already present at the codon level
  if (NROW(var_eitems[[target_type]]) > 0) {
    for (m in target_other) {
      if (NROW(var_eitems[[m]] > 0)) {
        assertable::assert_colnames(var_eitems[[m]], "GENOMIC_CHANGE",
                                    only_colnames = F, quiet = T)
        var_eitems[[m]] <- var_eitems[[m]] %>%
          dplyr::anti_join(
            dplyr::select(var_eitems[[target_type]], GENOMIC_CHANGE),
            by = c("GENOMIC_CHANGE"))
      }
    }
  }
  return(var_eitems)
}

log_var_eitem_stats <- function(var_eitems = NULL,
                               target_type = "exact"){

  invisible(
    assertthat::assert_that(!is.null(var_eitems),
                            msg = "Object 'var_eitems' cannot be NULL"))
  invisible(
    assertthat::assert_that(target_type == "exact" |
                              target_type == "codon" |
                              target_type == "exon" |
                              target_type == "gene",
                            msg = paste0("Argument 'target_type' can only",
                                         " take on values 'codon' or 'exact'",
                                         " or 'exon' or 'gene'")))


  rlogging::message(
    paste0("Found n = ",
           NROW(var_eitems[[target_type]]),
           " other clinical evidence item(s) at the ", target_type,
           " level, ",
           length(unique(var_eitems[[target_type]]$GENOMIC_CHANGE)),
           " unique variant(s)")
  )

  if (NROW(var_eitems[[target_type]]) > 0) {
    rlogging::message(
      paste0("Variants: ",
             paste(unique(paste(var_eitems[[target_type]]$SYMBOL,
                                var_eitems[[target_type]]$CONSEQUENCE,
                                var_eitems[[target_type]]$PROTEIN_CHANGE,
                                sep = ":")),
                   collapse = ", "))
    )
  }
}
