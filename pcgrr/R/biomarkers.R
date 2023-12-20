#' Function that matches clinical evidence items (CIVIC, CBMDB)
#' against somatic cancer variants detected in tumor
#'
#' @param sample_calls data frame with sample variants
#' @param annotation_tags list with annotation tags for display in report
#' @param eitems data frame with clinical evidence items
#'
#' @return list
#' @export

get_clin_assocs_snv_indel <- function(sample_calls,
                                      annotation_tags = NULL,
                                      eitems = NULL) {

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

    var_eitems[["exact"]] <- var_eitems[["exact"]] |>
      dplyr::bind_rows(var_eitems_exact) |>
      pcgrr::remove_cols_from_df(
        cnames = c("EITEM_CONSEQUENCE",
                   "EITEM_CODON",
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
    if (NROW(var_eitems_regional) > 0) {
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
  all_var_evidence_items <- all_var_evidence_items |>
    dplyr::bind_rows(var_eitems[["exact"]]) |>
    dplyr::bind_rows(var_eitems[["codon"]]) |>
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
  if (NROW(all_var_evidence_items) > 0) {
    variant_tags <-
      annotation_tags[["all"]][!annotation_tags[["all"]]
                                              %in% c("CIVIC_ID",
                                                     "CIVIC_ID_SEGMENT",
                                                     "CGI_ID_SEGMENT",
                                                     "CGI_ID")]
    variant_set <-  all_var_evidence_items |>
      dplyr::select(dplyr::one_of(variant_tags)) |>
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
#' @export
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

  assertable::assert_colnames(
    eitems,
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
      if (NROW(eitem_hits) > 0) {
        variant_set <- dplyr::bind_rows(variant_set, eitem_hits) |>
          dplyr::select(
            dplyr::one_of(annotation_tags[["cna_display"]])) |>
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

#' Function that loads specific set of clinical variant evidence items (CIViC + CGI) based
#' on given parameters (mutation type, variant origin, tumor type etc)
#'
#' @param eitems_raw complete set of clinical variant evidence items
#' @param ontology phenotype ontology data frame
#' @param alteration_types types of alteration ('MUT', 'CNA', 'MUT_LOF')
#' @param origin variant origin ('Somatic','Germline')
#' @param tumor_type_specificity tumor type specificity ('any', 'specific')
#' @param tumor_type primary tumor site
#'
#' @return eitems variant evidence items
#'
#'
#' @export
load_eitems <- function(eitems_raw = NULL,
                        ontology = NULL,
                        alteration_types = c("MUT"),
                        origin = "Somatic",
                        tumor_type_specificity = NULL,
                        tumor_type = NULL) {

  invisible(assertthat::assert_that(
    !is.null(eitems_raw),
    msg = "'eitems_raw' is NULL - existing"))

  invisible(assertthat::assert_that(
    !is.null(alteration_types),
    msg = "'alteration_types' is NULL - existing"))

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
      is.character(alteration_types),
      msg = "'alteration_types' must be a character vector, any combination of ('MUT','MUT_LOF','CNA')"
    )
  )



allowed_alt_types <- c("MUT", "CNA", "MUT_LOF")
assertthat::assert_that(
  all(alteration_types %in% allowed_alt_types),
  msg = paste0("Argument 'alteration_types' can only take the following values: ",
               paste0(allowed_alt_types, collapse = ", "), " and NOT: ",
               paste0(alteration_types[!alteration_types %in% allowed_alt_types], collapse = ", ")))


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
  eitems_all <- data.frame()

  for(alteration_type in alteration_types){
    eitems_alteration_type <-
      pcgrr::load_all_eitems(
        eitems_raw = eitems_raw,
        alteration_type = alteration_type,
        origin = origin)


    if (tumor_type_specificity == "any") {
      pcgrr::log4r_info(
        paste0(
          "Loading ", alteration_type, " biomarkers for precision oncology",
          " - any tumortype"))
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
        eitems_alteration_type, c("DISEASE_ONTOLOGY_ID"),
        only_colnames = F, quiet = T)
      assertable::assert_colnames(
        ontology,
        c("primary_site", "do_id", "cui", "cui_name"),
        only_colnames = F, quiet = T)

      eitems_alteration_type <-
        pcgrr::filter_eitems_by_site(
          eitems_alteration_type,
          ontology = ontology,
          primary_site = tumor_type)
      pcgrr::log4r_info(
        paste0("Loading ", alteration_type,
               " biomarkers for precision oncology - ",
               tumor_type))
    }

    eitems_all <- eitems_all |>
      dplyr::bind_rows(eitems_alteration_type)
  }
  return(eitems_all)


}

#' Function that loads all evidence items from CIViC and CGI, and
#' combines them in a unified data.frame
#'
#' @param eitems_raw raw data frame with evidence items
#' @param alteration_type type of alteration ('MUT','MUT_LOF','CNA')
#' @param origin variant origin ('Germline','Somatic')
#'
#' @return all_eitems
#'
#' @export
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
        eitems_raw[[db]] |>
          dplyr::filter(.data$ALTERATION_TYPE == alteration_type &
                          !is.na(.data$EITEM_CONSEQUENCE) &
                          stringr::str_detect(.data$VARIANT_ORIGIN, origin)) |>
          dplyr::rename(CNA_TYPE = .data$EITEM_CONSEQUENCE) |>
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
        eitems_raw[[db]] |>
        dplyr::filter(.data$ALTERATION_TYPE == alteration_type &
                        stringr::str_detect(.data$VARIANT_ORIGIN, origin)) |>
        pcgrr::remove_cols_from_df(
          cnames =
            c("VARIANT_NAME",
              "STATUS",
              "DRUG_INTERACTION_TYPE",
              "MOLECULE_CHEMBL_ID",
              "MAPPING_RANK",
              "BIOMARKER_DESCRIPTION"))
    }
    selected_eitems[[db]] <- selected_eitems[[db]] |>
      dplyr::mutate(SOURCE_DB = db)

  }

  all_eitems <- dplyr::bind_rows(selected_eitems[["civic"]],
                                 selected_eitems[["cgi"]])

  return(all_eitems)
}


#' Function that matches variants to evidence items
#'
#' param sample_calls data frame with variant calls
#' param db database with evidence items ('civic','cgi')
#' param colset character vector with column names to pull out from sample_calls
#' param eitems raw list of evidence items
#' param region_marker logical indication if region biomarkers are to be matched or not
#'
#'
#' export
# match_eitems_to_var <- function(sample_calls,
#                                db = "civic",
#                                colset = NULL,
#                                eitems = NULL,
#                                region_marker = T) {
#
#   invisible(assertthat::assert_that(
#     db == "civic" | db == "cgi",
#     msg = "Argument 'db' can be one of 'civic' or 'cgi'"))
#   invisible(assertthat::assert_that(
#     !is.null(eitems),
#     msg = "Argument 'eitems' cannot be NULL"))
#   invisible(assertthat::assert_that(
#     is.data.frame(eitems),
#     msg = paste0(
#       "Argument 'eitems' must be of type data frame, not ",
#       class(eitems))))
#   invisible(assertthat::assert_that(
#     !is.null(sample_calls) & is.data.frame(sample_calls),
#     msg = paste0(
#       "Argument 'sample_calls' must be of type data frame, not ",
#       class(sample_calls))))
#   assertable::assert_colnames(
#     eitems, c("EVIDENCE_ID", "SYMBOL","HGVS_ALIAS","SOURCE_DB"),
#     only_colnames = F, quiet = T)
#
#   invisible(assertthat::assert_that(!is.null(colset)))
#   invisible(assertthat::assert_that(is.character(colset)))
#
#   evidence_identifiers <- c("CIVIC_ID", "CIVIC_ID_SEGMENT")
#   if (region_marker == T) {
#     evidence_identifiers <- c("CIVIC_ID_SEGMENT", "CIVIC_ID")
#   }
#   eitems_db <- eitems |>
#     dplyr::filter(.data$SOURCE_DB == "civic") |>
#     dplyr::distinct()
#
#
#   if(db == "cgi"){
#     evidence_identifiers <- c("CGI_ID", "CGI_ID_SEGMENT")
#     if (region_marker == T) {
#       evidence_identifiers <- c("CGI_ID_SEGMENT", "CGI_ID")
#     }
#     eitems_db <- eitems |>
#       dplyr::filter(.data$SOURCE_DB == "cgi") |>
#       dplyr::distinct()
#   }
#
#   var_eitems <- list()
#   #var_eitems_exact <- data.frame()
#
#   assertable::assert_colnames(
#     sample_calls,
#     c(evidence_identifiers, colset),
#     only_colnames = F, quiet = T)
#
#   sample_calls_db <- sample_calls |>
#     dplyr::filter(!is.na(!!rlang::sym(evidence_identifiers[1])))
#   if (NROW(sample_calls_db) > 0) {
#     var_eitems_by_id <- as.data.frame(sample_calls_db |>
#       tidyr::separate_rows(!!rlang::sym(evidence_identifiers[1]), sep = ",") |>
#       dplyr::select(
#         dplyr::one_of(colset)) |>
#       dplyr::rename(EVIDENCE_ID = !!rlang::sym(evidence_identifiers[1])) |>
#       dplyr::mutate(EVIDENCE_ID = as.character(.data$EVIDENCE_ID)) |>
#       pcgrr::remove_cols_from_df(cnames = evidence_identifiers[2])
#     )
#
#     if (NROW(var_eitems_by_id) > 0) {
#       var_eitems[['by_id']] <- as.data.frame(
#         var_eitems_by_id |>
#           dplyr::inner_join(eitems_db,
#                            by = c("EVIDENCE_ID", "SYMBOL")) |>
#           dplyr::distinct() |>
#           pcgrr::remove_cols_from_df(
#             cnames = c("HGVS_ALIAS", evidence_identifiers))
#       )
#
#       if(NROW(var_eitems[['by_id']]) > 0){
#         var_eitems[['all']] <- var_eitems[['by_id']]
#       }
#
#     }
#   }
#
#   ## Add additional var_eitems based on matching against
#   ## HGVS (protein_change) + SYMBOL
#
#   if(region_marker == F){
#     eitems_hgvs <- eitems_db |>
#       dplyr::filter(!is.na(.data$HGVS_ALIAS))
#
#     if(NROW(eitems_hgvs) > 0){
#       eitems_hgvs <- eitems_hgvs |>
#         tidyr::separate_rows(
#           .data$HGVS_ALIAS, sep = "\\|") |>
#         dplyr::filter(
#           !stringr::str_detect(HGVS_ALIAS, "^rs")) |>
#         dplyr::rename(PROTEIN_CHANGE = .data$HGVS_ALIAS)
#
#       vars_hgvs_mapped <- sample_calls |>
#         dplyr::filter(!is.na(.data$PROTEIN_CHANGE)) |>
#         dplyr::select(dplyr::one_of(colset))
#
#       if(NROW(vars_hgvs_mapped) > 0){
#         var_eitems_hgvs_mapped <- as.data.frame(vars_hgvs_mapped |>
#           dplyr::inner_join(
#             eitems_hgvs, by = c("SYMBOL","PROTEIN_CHANGE")) |>
#           dplyr::distinct() |>
#           pcgrr::remove_cols_from_df(cnames = evidence_identifiers)
#         )
#
#         ## skip duplicate evidence items already found from
#         ## exact matching at genomic level
#         if(NROW(var_eitems_hgvs_mapped) > 0){
#           if(NROW(var_eitems[['by_id']]) > 0){
#             var_eitems_hgvs_mapped <-
#               var_eitems_hgvs_mapped |>
#               dplyr::anti_join(
#                 var_eitems[['by_id']], by = c("GENOMIC_CHANGE"))
#           }
#
#           if(NROW(var_eitems_hgvs_mapped) > 0){
#             var_eitems[['all']] <- var_eitems_exact |>
#               dplyr::bind_rows(var_eitems_hgvs_mapped) |>
#               dplyr::distinct()
#           }
#         }
#       }
#     }
#
#   }else {
#
#     ## Add additional var_eitems based on matching against
#     ## Refererence amino acid + position (e.g. codon) + SYMBOL
#
#     eitems_hgvs_codon <- eitems_db |>
#       dplyr::filter(!is.na(.data$HGVS_ALIAS)) |>
#       dplyr::filter(BIOMARKER_MAPPING == "codon")
#
#     if(NROW(eitems_hgvs_codon) > 0){
#       eitems_hgvs_codon <- eitems_hgvs_codon |>
#         tidyr::separate_rows(.data$HGVS_ALIAS, sep = "\\|") |>
#         dplyr::filter(
#           !stringr::str_detect(.data$HGVS_ALIAS, "^rs")) |>
#         dplyr::rename(AA_CODON = HGVS_ALIAS)
#
#       colset <- c('AA_CODON', colset)
#
#       vars_codon_mapped <- sample_calls |>
#         dplyr::filter(
#           !is.na(.data$PROTEIN_CHANGE) &
#             !is.na(AMINO_ACID_START) &
#             !is.na(AMINO_ACID_END) &
#             !is.na(Amino_acids) &
#             AMINO_ACID_START == AMINO_ACID_END) |>
#         dplyr::mutate(
#           AA_CODON = paste0(
#             stringr::str_replace(
#               Amino_acids, "/([A-Z]|\\*)$",""
#             ), AMINO_ACID_START
#           )) |>
#         dplyr::select(dplyr::one_of(colset))
#
#       if(NROW(vars_codon_mapped) > 0){
#         var_eitems_codon_mapped <- as.data.frame(
#           vars_codon_mapped |>
#             dplyr::inner_join(
#               eitems_hgvs_codon, by = c("SYMBOL","AA_CODON")) |>
#             dplyr::distinct() |>
#             pcgrr::remove_cols_from_df(
#               cnames = evidence_identifiers)
#         )
#
#         ## skip duplicate evidence items already found from
#         ## exact matching at genomic level
#         if(nrow(var_eitems_codon_mapped) > 0){
#           if(NROW(var_eitems[['by_id']]) > 0){
#             var_eitems_codon_mapped <- var_eitems_codon_mapped |>
#               dplyr::select(-c("AA_CODON")) |>
#               dplyr::anti_join(
#                 var_eitems[['by_id']],
#                 by = c("GENOMIC_CHANGE","BIOMARKER_MAPPING"))
#           }
#
#           var_eitems[['all']] <- var_eitems[['all']] |>
#             dplyr::bind_rows(var_eitems_codon_mapped) |>
#             dplyr::distinct()
#         }
#       }
#     }
#
#   }
#
#   return(var_eitems[['all']])
#
# }

#' Function that makes a quality control check of evidence items assigned
#' to variants
#'
#' @param var_eitems variant-evidence items
#' @param marker_type type of biomarker
#'
#' @export
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
      filtered_var_eitems <- var_eitems |>
        dplyr::filter(!is.na(.data$EITEM_CODON) & .data$BIOMARKER_MAPPING == "codon") |>
        dplyr::filter(.data$EITEM_CODON <= .data$AMINO_ACID_END &
                        .data$EITEM_CODON >= .data$AMINO_ACID_START &
                        (!is.na(.data$EITEM_CONSEQUENCE) &
                           startsWith(.data$CONSEQUENCE, .data$EITEM_CONSEQUENCE) |
                           is.na(.data$EITEM_CONSEQUENCE)) &
                        .data$CODING_STATUS == "coding")

    }
  }

  if (marker_type == "exon") {
    if (nrow(var_eitems[!is.na(var_eitems$EITEM_EXON) &
                    var_eitems$BIOMARKER_MAPPING == "exon", ]) > 0) {
      filtered_var_eitems <- var_eitems |>
        dplyr::filter(!is.na(.data$EITEM_EXON) & .data$BIOMARKER_MAPPING == "exon") |>
        dplyr::filter(.data$AFFECTED_EXON == .data$EITEM_EXON &
                        (!is.na(.data$EITEM_CONSEQUENCE) &
                          startsWith(.data$CONSEQUENCE, .data$EITEM_CONSEQUENCE) |
                         is.na(.data$EITEM_CONSEQUENCE)) & .data$CODING_STATUS == "coding")
    }
  }

  if (marker_type == "gene") {
    if (nrow(var_eitems[var_eitems$BIOMARKER_MAPPING == "gene", ]) > 0) {
      filtered_var_eitems <- var_eitems |>
        dplyr::filter(.data$BIOMARKER_MAPPING == "gene" &
                        .data$CODING_STATUS == "coding") |>
        dplyr::filter((!is.na(.data$EITEM_CONSEQUENCE) &
                         startsWith(.data$CONSEQUENCE, .data$EITEM_CONSEQUENCE) |
                         is.na(.data$EITEM_CONSEQUENCE)) & .data$CODING_STATUS == "coding")

    }
  }

  if (nrow(filtered_var_eitems) > 0) {

    if("LOSS_OF_FUNCTION" %in% colnames(filtered_var_eitems) &
       "ALTERATION_TYPE" %in% colnames(filtered_var_eitems)){

      filtered_var_eitems <- filtered_var_eitems |>
        dplyr::filter((.data$LOSS_OF_FUNCTION == T &
                         .data$ALTERATION_TYPE == "MUT_LOF") |
                        is.na(.data$LOSS_OF_FUNCTION) |
                        (.data$LOSS_OF_FUNCTION == F &
                           .data$ALTERATION_TYPE != "MUT_LOF"))
    }
  }

  filtered_var_eitems <- filtered_var_eitems |>
    pcgrr::remove_cols_from_df(cnames = c("EITEM_CONSEQUENCE",
                                          "EITEM_CODON",
                                          "EITEM_EXON"))

  return(filtered_var_eitems)

}

#' Function that filters clinical evidence items by tumor type/primary site
#'
#' @param eitems data frame with clinical evidence items
#' @param ontology phenotype ontology data frame
#' @param primary_site primary tumor site
#'
#' @export
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
                    .data$primary_site, .data$do_id, .data$cui, .data$cui_name),
      data.frame("primary_site" = primary_site, stringsAsFactors = F),
      by = c("primary_site" = "primary_site")) |>
    dplyr::filter(!is.na(.data$do_id)) |>
    dplyr::distinct()

  eitems <- eitems |>
    dplyr::semi_join(tumor_phenotypes_site,
                     by = c("DISEASE_ONTOLOGY_ID" = "do_id"))

  return(eitems)
}

#' Function that structures variant evidence items according
#' to strength of evidence
#'
#' @param var_eitems variant evidence items
#' @param annotation_tags annotation tags to include for display
#' @param alteration_type type of alteration ('MUT','CNA')
#'

#' @export
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
      clin_eitems_list[[type]][["any"]] <- var_eitems |>
        dplyr::select(dplyr::one_of(tags_display)) |>
        dplyr::filter(.data$EVIDENCE_TYPE == stringr::str_to_title(type)) |>
        dplyr::arrange(.data$EVIDENCE_LEVEL, dplyr::desc(.data$RATING))
      if (nrow(clin_eitems_list[[type]][["any"]]) > 0) {
        clin_eitems_list[[type]][["A_B"]] <-
          clin_eitems_list[[type]][["any"]] |>
          dplyr::filter(stringr::str_detect(.data$EVIDENCE_LEVEL, "^(A|B|B1|B2):"))

        if (NROW(clin_eitems_list[[type]][["A_B"]]) > 0) {
          clin_eitems_list[[type]][["A_B"]] <-
            clin_eitems_list[[type]][["A_B"]] |>
            dplyr::arrange(.data$EVIDENCE_LEVEL, dplyr::desc(.data$RATING))
        }

        clin_eitems_list[[type]][["C_D_E"]] <-
          clin_eitems_list[[type]][["any"]] |>
          dplyr::filter(stringr::str_detect(.data$EVIDENCE_LEVEL, "^(C|D|E):"))

        if (NROW(clin_eitems_list[[type]][["C_D_E"]]) > 0) {
          clin_eitems_list[[type]][["C_D_E"]] <-
            clin_eitems_list[[type]][["C_D_E"]] |>
            dplyr::arrange(.data$EVIDENCE_LEVEL, dplyr::desc(.data$RATING))
        }
      }
    }
  }
  return(clin_eitems_list)


}


#' Function that removes redundancy in variant evidence items (i.e. if
#' a variant is assicated with evidence at the codon level, evidence
#' at the exon/gene level is ignored)
#'
#' @param var_eitems data frame with variant evidence items
#' @param target_type which resolution level should be used as the
#' "best" level ('exact' or 'codon)
#' @param target_other resolution levels for other evidence items
#' that should be ignored if evidence is found at the target_type level
#'
#'
#' @export
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
    assertable::assert_colnames(
      var_eitems[[target_type]], "GENOMIC_CHANGE",
      only_colnames = F, quiet = T)
    for (m in target_other) {
      if (NROW(var_eitems[[m]]) > 0) {
        assertable::assert_colnames(var_eitems[[m]], "GENOMIC_CHANGE",
                                    only_colnames = F, quiet = T)
        var_eitems[[m]] <- var_eitems[[m]] |>
          dplyr::anti_join(
            dplyr::select(var_eitems[[target_type]], .data$GENOMIC_CHANGE),
            by = c("GENOMIC_CHANGE"))
      }
    }
  }
  return(var_eitems)
}

#' Function that logs the number of evidence items found, for different
#' levels of resolution
#'
#' @param var_eitems data frame with variant evidence items
#' @param target_type resolution of evidence items
#'
#'
#'

#' @export
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


  pcgrr::log4r_info(
    paste0("Found n = ",
           NROW(var_eitems[[target_type]]),
           " clinical evidence item(s) at the ", target_type,
           " level, ",
           length(unique(var_eitems[[target_type]]$GENOMIC_CHANGE)),
           " unique variant(s)")
  )

  if (NROW(var_eitems[[target_type]]) > 0) {
    assertable::assert_colnames(
      var_eitems[[target_type]],
      c("SYMBOL","CONSEQUENCE","PROTEIN_CHANGE"),
      only_colnames = F, quiet = T)

    variants_found_log <-
      paste(unique(paste(var_eitems[[target_type]]$SYMBOL,
                         var_eitems[[target_type]]$CONSEQUENCE,
                         var_eitems[[target_type]]$PROTEIN_CHANGE,
                         sep = ":")),
            collapse = ", ")
    if(nchar(variants_found_log) <= 200){
      pcgrr::log4r_info(
        variants_found_log
      )
    }
  }
}

#' Function that expands biomarker evidence items with variant annotations
#'
#' @param callset list object with 'variant' and 'biomarker_evidence' data
#' frames
#' @param variant_origin 'somatic' or 'germline'
#' @param target_genes data frame with target genes of interest
#'
#' @export
#'
expand_biomarker_items <- function(
    callset = NULL,
    variant_origin = "somatic",
    target_genes = NULL){

  if("variant" %in% names(callset) &
     "biomarker_evidence" %in% names(callset)){

    variant_properties <-
      c("VAR_ID",
        "GENOMIC_CHANGE",
        "GENOME_VERSION",
        "SAMPLE_ID",
        "GENOTYPE",
        "VARIANT_CLASS",
        "SYMBOL",
        "GENENAME",
        "ENTREZGENE",
        "CONSEQUENCE",
        "PROTEIN_CHANGE",
        "MUTATION_HOTSPOT",
        "CDS_CHANGE",
        "LOSS_OF_FUNCTION",
        "HGVSc",
        "HGVSp",
        "REFSEQ",
        "OFFICIAL_GENENAME",
        "PREDICTED_EFFECT",
        "PROTEIN_DOMAIN",
        "DBSNP",
        "CLINVAR",
        "COSMIC",
        "VEP_ALL_CSQ")

    if(variant_origin == "germline"){
      variant_properties <- c(
        variant_properties,
        "CLINVAR_CLASSIFICATION",
        "CPSR_CLASSIFICATION"
      )
    }
    if(variant_origin == "somatic"){
      variant_properties <- c(
        variant_properties,
        "DP_TUMOR",
        "AF_TUMOR",
        "DP_CONTROL",
        "AF_CONTROL"
      )
    }

    ## check col existence callset[['variant]], variant_properties

    for (type in c(pcgrr::evidence_types,
                   "all")) {
      for (elevel in c("any", "A_B", "C_D_E")) {
        if(NROW(callset[['biomarker_evidence']][[type]][[elevel]]) > 0){
          callset[['biomarker_evidence']][[type]][[elevel]] <-
            callset[['biomarker_evidence']][[type]][[elevel]] |>
            dplyr::left_join(
              dplyr::select(
                callset[['variant']],
                variant_properties),
              by = c("VAR_ID")) |>
            dplyr::arrange(
              .data$EVIDENCE_LEVEL,
              .data$PROTEIN_CHANGE,
              dplyr::desc(
                .data$RATING))

          if(variant_origin == "germline"){
            callset[['biomarker_evidence']][[type]][[elevel]] <-
              callset[['biomarker_evidence']][[type]][[elevel]] |>
              dplyr::filter(
                (!is.na(CLINVAR_CLASSIFICATION) &
                   stringr::str_detect(
                     tolower(CLINVAR_CLASSIFICATION), "pathogenic")) |
                  (is.na(CLINVAR_CLASSIFICATION) &
                     !is.na(CPSR_CLASSIFICATION) &
                     stringr::str_detect(
                       tolower(CPSR_CLASSIFICATION), "pathogenic"))
              )

            if(NROW(callset[['biomarker_evidence']][[type]][[elevel]]) > 0 &
               is.data.frame(target_genes) &
               NROW(target_genes) > 0 &
               "ENTREZGENE" %in% colnames(target_genes)){
              callset[['biomarker_evidence']][[type]][[elevel]] <-
                callset[['biomarker_evidence']][[type]][[elevel]] |>
                dplyr::semi_join(target_genes, by = "ENTREZGENE")

            }
          }
        }
      }
    }
  }

  return(callset)

}
