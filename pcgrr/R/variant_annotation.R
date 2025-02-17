#' Function that appends multiple HTML annotation links to variant identifiers
#' e.g. COSMIC, CLINVAR, REFSEQ etc
#'
#' @param var_data_df data frame with variant entries
#' @param vartype 'snv_indel' or 'cna', or 'exp'
#' @param skip elements to be ignored during annotation
#'
#' @export
append_annotation_links <- function(
    var_data_df,
    vartype = "snv_indel",
    skip = NULL) {

  i <- 1
  while (i <= nrow(pcgrr::variant_db_url)) {

    name <-
      pcgrr::variant_db_url[i, ]$name
    link_display_var <-
      pcgrr::variant_db_url[i, ]$link_display_var
    rename <-
      pcgrr::variant_db_url[i, ]$rename

    if(rename == T &
       name %in% colnames(var_data_df) &
       !(link_display_var %in% colnames(var_data_df))){
      names(var_data_df)[names(var_data_df) == name] <-
        link_display_var
    }

    if (!name %in% skip) {
      group_by_var <- c("VAR_ID","ENTREZGENE")
      if(vartype == "cna"){
        group_by_var <- c("VAR_ID",
                          "ENTREZGENE",
                          "TRANSCRIPT_START",
                          "TRANSCRIPT_END")
      }
      if(vartype == "exp"){
        group_by_var <- "ENTREZGENE"
      }
      url_prefix <- pcgrr::variant_db_url[i, ]$url_prefix
      link_key_var <- pcgrr::variant_db_url[i, ]$link_key_var
      link_display_var <- pcgrr::variant_db_url[i, ]$link_display_var
      if (!(name %in% colnames(var_data_df)) &
          link_key_var %in% colnames(var_data_df) &
          link_display_var %in% colnames(var_data_df) &
          length(unique(group_by_var %in% colnames(var_data_df))) == 1) {
        annotation_links <-
          pcgrr::generate_annotation_link(
            var_df = var_data_df,
            vardb = name,
            vartype = vartype,
            group_by_var = group_by_var,
            url_prefix = url_prefix,
            link_key_var = link_key_var,
            link_display_var = link_display_var
          )
        if (nrow(annotation_links) > 0) {
          if(vartype == "cna"){
            var_data_df <- var_data_df |>
              dplyr::left_join(
                dplyr::rename(
                  annotation_links,
                  !!rlang::sym(name) := "link"),
                by = c("VAR_ID","ENTREZGENE",
                       "TRANSCRIPT_START",
                       "TRANSCRIPT_END"))
          }else{
            if(vartype == "exp"){
              var_data_df <- var_data_df |>
                dplyr::left_join(
                  dplyr::rename(
                    annotation_links,
                    !!rlang::sym(name) := "link"),
                  by = c("ENTREZGENE"))
            }
            if(vartype == "snv_indel"){
              var_data_df <- var_data_df |>
                dplyr::left_join(
                  dplyr::rename(
                    annotation_links,
                    !!rlang::sym(name) := "link"),
                  by = c("VAR_ID","ENTREZGENE"))
            }
          }
        }else{
          var_data_df[, name] <- NA
        }
      }
    }
    i <- i + 1
  }
  return(var_data_df)
}


#' Function that adds GWAS citation/phenotype to GWAS hit found through PCGR annotation
#'
#' @param vcf_data_df Data frame of sample variants from VCF
#' @param ref_data PCGR/CPSR reference data object
#'
#' @return vcf_data_df
#'
#' @export
append_gwas_citation_phenotype <-
  function(vcf_data_df = NULL,
           ref_data = NULL){


    invisible(assertthat::assert_that(
      !is.null(vcf_data_df),
      msg = "Argument 'vcf_data_df' cannot not be NULL"))
    invisible(assertthat::assert_that(
      is.data.frame(vcf_data_df),
      msg = paste0("Argument 'vcf_data_df' must be of type 'data.frame'")))
    invisible(assertthat::assert_that(
      !is.null(ref_data),
      msg = "Argument 'ref_data' cannot not be NULL"))
    invisible(assertthat::assert_that(
      !is.null(ref_data$variant),
      msg = "'ref_data$variant' cannot not be NULL"))
    invisible(assertthat::assert_that(
      is.data.frame(ref_data$variant$gwas),
      msg =
        "Argument 'ref_data$variant$gwas' must be of type 'data.frame'"))
    assertable::assert_colnames(
      vcf_data_df, c("GWAS_HIT", "VAR_ID"),
      only_colnames = F, quiet = T)
    assertable::assert_colnames(
      ref_data$variant$gwas,
      c("GWAS_HIT",
        "CITATION",
        "GWAS_CITATION_EXP",
        "GWAS_PHENOTYPE"),
      only_colnames = F, quiet = T)

    if ("GWAS_HIT" %in% colnames(vcf_data_df) &
        "VAR_ID" %in% colnames(vcf_data_df)) {

      pcgrr::log4r_info(paste0("Adding citations/phenotypes underlying ",
                               "GWAS hits (NHGRI-EBI GWAS Catalog)"))
      feature_df <- dplyr::select(
        vcf_data_df,
        c("GWAS_HIT","VAR_ID")) |>
        dplyr::filter(
          !is.na(.data$GWAS_HIT)) |>
        dplyr::distinct()
      if (nrow(feature_df) == 0) {
        vcf_data_df$GWAS_CITATION <- NA
        vcf_data_df$GWAS_CITATION2 <- NA
        vcf_data_df$GWAS_PHENOTYPE <- NA
        return(vcf_data_df)
      }

      feature_df <- as.data.frame(
        feature_df |>
          tidyr::separate_rows(.data$GWAS_HIT, sep = ",") |>
          dplyr::left_join(
            dplyr::select(
              ref_data$variant$gwas,
              c("GWAS_HIT",
                "GWAS_CITATION_EXP",
                "CITATION",
                "GWAS_PHENOTYPE")),
            by = "GWAS_HIT") |>
          dplyr::filter(
            !is.na(.data$GWAS_CITATION_EXP))
      )

      if (nrow(feature_df) == 0) {
        vcf_data_df$GWAS_CITATION <- NA
        vcf_data_df$GWAS_CITATION2 <- NA
        vcf_data_df$GWAS_PHENOTYPE <- NA
        return(vcf_data_df)
      }

      feature_df <- as.data.frame(
        feature_df |>
          dplyr::group_by(.data$VAR_ID) |>
          dplyr::summarise(
            GWAS_PHENOTYPE = paste(
              unique(.data$GWAS_PHENOTYPE), collapse = "; "),
            GWAS_CITATION = paste(
              unique(.data$GWAS_CITATION_EXP), collapse = "; "),
            GWAS_CITATION2 = paste(
              unique(.data$GWAS_CITATION), collapse = "; ")) |>
          dplyr::mutate(
            GWAS_CITATION2 = dplyr::if_else(
              .data$GWAS_CITATION2 == "NA",
              as.character(NA),
              as.character(.data$GWAS_CITATION2))) |>
          dplyr::mutate(
            GWAS_CITATION = dplyr::if_else(
              .data$GWAS_CITATION == "NA",
              as.character(NA),
              as.character(.data$GWAS_CITATION))) |>
          dplyr::mutate(
            GWAS_PHENOTYPE = dplyr::if_else(
              .data$GWAS_PHENOTYPE == "NA",
              as.character(NA),
              as.character(.data$GWAS_PHENOTYPE))) |>
          dplyr::filter(
            !is.na(.data$GWAS_CITATION) &
              !is.na(.data$GWAS_PHENOTYPE))
      )
      if (nrow(feature_df) > 0) {

        pcgrr::log4r_info(paste0(
          "Found n = ",
          nrow(feature_df),
          " variants associated with genome-wide association studies"))

        vcf_data_df <- vcf_data_df |>
          dplyr::left_join(feature_df, by = c("VAR_ID" = "VAR_ID"))
      }else{
        vcf_data_df$GWAS_CITATION <- NA
        vcf_data_df$GWAS_CITATION2 <- NA
        vcf_data_df$GWAS_PHENOTYPE <- NA
      }
    }
  }


#' Function that adds TCGA annotations (cohort, frequency etc.) to variant identifiers
#'
#' @param var_df data frame with variants
#' @param linktype type of link
#'
#' @return var_df
#'
#' @export
append_tcga_var_link <- function(var_df,
                                 linktype = "dbsource") {


  if (any(grepl(paste0("^TCGA_FREQUENCY$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df)))) {

    var_df_unique_slim <- dplyr::select(
      var_df, .data$VAR_ID, .data$TCGA_FREQUENCY) |>
      dplyr::filter(!is.na(.data$TCGA_FREQUENCY)) |>
      dplyr::distinct()
    if (NROW(var_df_unique_slim) > 0) {
      var_df_unique_slim_melted <- var_df_unique_slim |>
        tidyr::separate_rows(.data$TCGA_FREQUENCY, sep = ",") |>
        tidyr::separate(
          .data$TCGA_FREQUENCY,
          c("tcga_cancer_code", "percentage",
            "affected", "cohort_size"),
          sep = "\\|", convert = T) |>
        dplyr::left_join(
          pcgrr::tcga_cohorts, by = "tcga_cancer_code") |>
        dplyr::arrange(
          .data$VAR_ID,
          dplyr::desc(.data$percentage))
      if (linktype == "dbsource") {
        var_df_unique_slim_melted <- var_df_unique_slim_melted |>
          dplyr::mutate(tmp_assoc = paste0(
            "<a href='https://portal.gdc.cancer.gov/projects/TCGA-",
            .data$tcga_cancer_code, "' target=\"_blank\">",
            .data$tcga_cancer_name, "</a>: ",
            .data$percentage, "% (",
            .data$affected, "/", .data$cohort_size, ")"))
      }

      var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) |>
        dplyr::summarise(TCGALINK = unlist(paste(.data$tmp_assoc, collapse = ", ")),
                         .groups = "drop") |>
        dplyr::select(.data$VAR_ID, .data$TCGALINK) |>
        dplyr::rename(TCGA_FREQUENCY = "TCGALINK")
        #magrittr::set_colnames(c("VAR_ID", "TCGA_FREQUENCY"))
      var_df <- dplyr::rename(
        var_df, TCGA_FREQUENCY_RAW = .data$TCGA_FREQUENCY)
      var_df <- dplyr::left_join(var_df, var_df_links,
                                 by = c("VAR_ID" = "VAR_ID"))
    }else{
      var_df$TCGA_FREQUENCY_RAW <- var_df$TCGA_FREQUENCY
    }
  }
  return(var_df)
}

#' Function that adds TFBS annotations (dbMTS) to genetic variant identifiers
#'
#' @param var_df data frame with variants
#'
#' @export
append_tfbs_annotation <-
  function(var_df) {

    if (any(grepl(paste0("^CONSEQUENCE$"), names(var_df))) &
        any(grepl(paste0("^VAR_ID$"), names(var_df))) &
        any(grepl(paste0("^REGULATORY_ANNOTATION$"), names(var_df)))) {

      var_df_unique_slim <-
        dplyr::select(var_df, .data$VAR_ID,
                      .data$REGULATORY_ANNOTATION,
                      .data$CONSEQUENCE) |>
        dplyr::filter(!is.na(.data$REGULATORY_ANNOTATION) &
                        stringr::str_detect(
                          .data$CONSEQUENCE,
                          "5_prime|upstream"
                        )) |>
        dplyr::distinct()

      if (nrow(var_df_unique_slim) > 0) {
        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim |>
            tidyr::separate_rows(.data$REGULATORY_ANNOTATION, sep=",") |>
            dplyr::filter(
              stringr::str_detect(
                .data$REGULATORY_ANNOTATION, "TF_binding_site_variant"
              ))
        )

        if (nrow(var_df_unique_slim_melted) > 0) {

          pcgrr::log4r_info(paste0(
            "Found TF binding site annotations for ",
            nrow(var_df_unique_slim)," variants"))

          var_df_unique_slim_melted <- as.data.frame(
            var_df_unique_slim_melted |>
              dplyr::mutate(
                REGULATORY_ANNOTATION = stringr::str_replace(
                  .data$REGULATORY_ANNOTATION,
                  "TF_binding_site_variant\\|MotifFeature\\|ENSM0[0-9]{1,}\\|",
                  "")
              ) |>
              tidyr::separate(.data$REGULATORY_ANNOTATION,
                              into = c('cons','matrix','motif_pos',
                                       'high_inf_pos','motif_score_change',
                                       'transcription_factors'),
                              sep = "\\|",
                              remove = T) |>

              tidyr::separate_rows(.data$transcription_factors) |>
              dplyr::mutate(
                TF_BINDING_SITE_VARIANT = dplyr::case_when(
                  .data$high_inf_pos == "N" ~ "Overlap: non-critical motif position",
                  .data$high_inf_pos == "Y" ~ "Overlap: critical motif position",
                  TRUE ~ as.character(NA)
                )
              ) |>
              dplyr::mutate(
                TF_BINDING_SITE_VARIANT_INFO =
                  paste(.data$transcription_factors, .data$matrix,
                        .data$motif_pos, .data$motif_score_change,
                        .data$high_inf_pos, sep="|")
              )
          )

          var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) |>
            dplyr::summarise(
              TF_BINDING_SITE_VARIANT = paste(unique(sort(.data$TF_BINDING_SITE_VARIANT)),
                                              collapse = ", "),
              TF_BINDING_SITE_VARIANT_INFO = paste(unique(
                .data$TF_BINDING_SITE_VARIANT_INFO),
                collapse = ", "),
              .groups = "drop") |>
            dplyr::select(.data$VAR_ID, .data$TF_BINDING_SITE_VARIANT,
                          .data$TF_BINDING_SITE_VARIANT_INFO)

          var_df <- dplyr::left_join(var_df, var_df_links,
                                     by = c("VAR_ID" = "VAR_ID"))

        }else{
          var_df$TF_BINDING_SITE_VARIANT <- NA
          var_df$TF_BINDING_SITE_VARIANT_INFO <- NA
        }
      }else{
        var_df$TF_BINDING_SITE_VARIANT <- NA
        var_df$TF_BINDING_SITE_VARIANT_INFO <- NA
      }
    }
    return(var_df)
  }

#' Function that adds miRNA target annotations (dbMTS) to
#' genetic variant identifiers
#'
#' @param var_df data frame with variants
#'
#' @export
append_dbmts_var_link <-
  function(var_df) {


    if (any(grepl(paste0("^DBMTS$"), names(var_df))) &
        any(grepl(paste0("^VAR_ID$"), names(var_df))) &
        any(grepl(paste0("^ENSEMBL_TRANSCRIPT_ID$"), names(var_df)))) {

      var_df_unique_slim <-
        dplyr::select(
          var_df, .data$VAR_ID, .data$CLINVAR_CLASSIFICATION,
          .data$DBMTS, .data$ENSEMBL_TRANSCRIPT_ID) |>
        dplyr::filter(!is.na(.data$DBMTS) &
                        !is.na(.data$ENSEMBL_TRANSCRIPT_ID)) |>
        dplyr::distinct()
      if (nrow(var_df_unique_slim) > 0) {

        var_df_unique_slim_melted <- as.data.frame(
          var_df_unique_slim |>
            tidyr::separate_rows(.data$DBMTS, sep = ",") |>
            tidyr::separate(.data$DBMTS, c("ens_trans_id", "mirna_id",
                                           "algorithms", "algorithms_call",
                                           "consensus_call"),
                            sep = "\\|", convert = T) |>
            dplyr::filter(.data$ens_trans_id == .data$ENSEMBL_TRANSCRIPT_ID)
        )
        if (nrow(var_df_unique_slim_melted) > 0) {
          var_df_unique_slim_melted <- var_df_unique_slim_melted |>
            dplyr::select(-c(.data$ENSEMBL_TRANSCRIPT_ID, .data$algorithms_call)) |>
            dplyr::mutate(miRNA_TARGET_HIT = dplyr::case_when(
              .data$consensus_call == "G" ~ "gain",
              .data$consensus_call == "L" ~ "loss",
              TRUE ~ as.character(NA)
            )) |>
            dplyr::mutate(
              algorithms = stringr::str_replace_all(
                stringr::str_replace(
                  stringr::str_replace(
                    stringr::str_replace(
                      .data$algorithms, "R","RNAHybrid"),
                    "TS","TargetScan"),
                  "M","miRanda"),
                "&"," / ")
            ) |>
            dplyr::mutate(
              miRNA_TARGET_HIT_PREDICTION =
                paste0("<a href='http://www.mirbase.org/cgi-bin/mirna_entry.pl?id",
                       "=",.data$mirna_id,"' target='_blank'>",.data$mirna_id,"</a> - ", .data$miRNA_TARGET_HIT,
                       " (",.data$algorithms,")")
            )


          var_df_links <- dplyr::group_by(var_df_unique_slim_melted, .data$VAR_ID) |>
            dplyr::summarise(
              miRNA_TARGET_HIT_PREDICTION = paste(.data$miRNA_TARGET_HIT_PREDICTION,
                                                  collapse = ", "),
              miRNA_TARGET_HIT = paste(unique(.data$miRNA_TARGET_HIT),
                                       collapse = ", "),
              .groups = "drop") |>
            dplyr::select(.data$VAR_ID, .data$miRNA_TARGET_HIT, .data$miRNA_TARGET_HIT_PREDICTION)

          var_df <- dplyr::left_join(var_df, var_df_links,
                                     by = c("VAR_ID" = "VAR_ID"))
        }else{
          var_df$miRNA_TARGET_HIT_PREDICTION <- NA
          var_df$miRNA_TARGET_HIT <- NA
        }
      }else{
        var_df$miRNA_TARGET_HIT_PREDICTION <- NA
        var_df$miRNA_TARGET_HIT <- NA
      }
    }

    return(var_df)
  }

#' Function that assigns HTML links to dbNSFP prediction entries
#'
#' @param var_df data frame with variant entries
#'
#' @return var_df
#'
#' @export
append_dbnsfp_var_link <- function(var_df) {

  if (any(grepl(paste0("EFFECT_PREDICTIONS"), names(var_df)))) {
    var_df <- var_df |>
      dplyr::mutate(PREDICTED_EFFECT = .data$EFFECT_PREDICTIONS) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":D,", ":Damaging,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":T,", ":Tolerated,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":SN,", ":SplicingNeutral,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":SN", ":SplicingNeutral"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":AS,", ":AffectSplicing,"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":AS", ":AffectSplicing"
      )) |>
      dplyr::mutate(PREDICTED_EFFECT = stringr::str_replace_all(
        .data$PREDICTED_EFFECT, ":PD,", ":ProbablyDamaging,"
      ))
    i <- 1
    while (i <= nrow(pcgrr::effect_prediction_algos)) {
      str_to_replace <-
        paste0(pcgrr::effect_prediction_algos[i, "algorithm"], ":")
      replacement_str <-
        paste0("<a href='",
               pcgrr::effect_prediction_algos[i, "url"], "' target='_blank'>",
               pcgrr::effect_prediction_algos[i, "display_name"], "</a>:")
      algorithm_display <-
        paste0(pcgrr::effect_prediction_algos[i, "display_name"], ":")
      var_df <- var_df |>
        dplyr::mutate(
          PREDICTED_EFFECT =
            stringr::str_replace(.data$PREDICTED_EFFECT,
                                 str_to_replace, replacement_str))
      i <- i + 1
    }
  }
  # else{
  #   var_df$PREDICTED_EFFECT <- NA
  # }
  return(var_df)

}

#' Function that adds link to targeted drugs (on and off-label)
#' for a list of variants and associated targeted
#'
#' @param var_df data frame with variants
#' @param primary_site primary site of tumor sample
#' @param ref_data PCGR/CPSR reference data object
#'
#' @export
append_targeted_drug_annotations <- function(
    var_df,
    primary_site = "Breast",
    ref_data = NULL) {

  if (any(grepl(paste0("^SYMBOL$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df))) &
      !is.null(ref_data)) {
    var_drug_df <- dplyr::select(
      var_df, c("VAR_ID","SYMBOL")) |>
      dplyr::filter(!is.na(.data$SYMBOL)) |>
      dplyr::distinct()
    if (nrow(var_drug_df) > 0) {
      var_drug_df <- var_drug_df |>
        dplyr::left_join(
          ref_data[['drug']][['inhibitors_any_label']],
          by = c("SYMBOL")) |>
        dplyr::left_join(
          dplyr::filter(
            ref_data[['drug']][['inhibitors_on_label']],
            .data$QUERY_SITE == primary_site),
          by = c("SYMBOL")) |>
        dplyr::select(-c("QUERY_SITE")) |>
        dplyr::distinct() |>
        dplyr::filter(
          !is.na(.data$TARGETED_INHIBITORS_ALL2)) |>
        dplyr::distinct()
      if (NROW(var_drug_df) > 0) {
        var_df <-
          dplyr::left_join(
            var_df,
            var_drug_df,
            by = c("VAR_ID","SYMBOL")
          ) |>
          dplyr::select(
            -c("TARGETED_INHIBITORS",
               "TARGETED_INHIBITORS_ALL")
          )
      }else{
        #vcf_data_df$TARGETED_INHIBITORS_OFF_LABEL <- NA
        var_df$TARGETED_INHIBITORS2 <- NA
        var_df$TARGETED_INHIBITORS_ALL2 <- NA
      }
    }
    else{
      var_df$TARGETED_INHIBITORS2 <- NA
      var_df$TARGETED_INHIBITORS_ALL2 <- NA
    }
  }

  return(var_df)
}

#' Function that adds link to targeted drugs (on and off-label)
#' for a list of variants and associated targeted
#'
#' @param var_df data frame with variants
#' @param primary_site primary tumor sample site
#' @param ref_data PCGR/CPSR reference data object
#'
#' @export
append_drug_var_link <- function(
    var_df,
    primary_site = "Breast",
    ref_data = NULL) {

  if (any(grepl(paste0("^SYMBOL$"), names(var_df))) &
      any(grepl(paste0("^VAR_ID$"), names(var_df))) &
      !is.null(ref_data)) {
    var_drug_df <- dplyr::select(
      var_df, c("VAR_ID","SYMBOL")) |>
      dplyr::filter(!is.na(.data$SYMBOL)) |>
      dplyr::distinct()
    if (nrow(var_drug_df) > 0) {
      var_drug_df <- var_drug_df |>
        dplyr::left_join(
            ref_data[['drug']][['inhibitors_any_label']],
          by = c("SYMBOL")) |>
        dplyr::left_join(
          dplyr::filter(
            ref_data[['drug']][['inhibitors_on_label']],
            .data$QUERY_SITE == primary_site),
          by = c("SYMBOL")) |>
        dplyr::select(-c("QUERY_SITE")) |>
        dplyr::distinct() |>
        dplyr::filter(
          !is.na(.data$TARGETED_INHIBITORS_ALL2)) |>
        dplyr::distinct()
      if (NROW(var_drug_df) > 0) {
        var_df <-
          dplyr::left_join(
            var_df,
            var_drug_df,
            by = c("VAR_ID","SYMBOL")
          )
      }else{
        #vcf_data_df$TARGETED_INHIBITORS_OFF_LABEL <- NA
        var_df$TARGETED_INHIBITORS <- NA
        var_df$TARGETED_INHIBITORS_ALL <- NA
      }
    }
    else{
      var_df$TARGETED_INHIBITORS <- NA
      var_df$TARGETED_INHIBITORS_ALL <- NA
    }
  }

  return(var_df)
}

#' Function that appends cancer gene evidence links
#'
#' @param var_data_df Data frame of sample variants from VCF
#' @param ref_data PCGR/CPSR reference data object
#' @param primary_site Primary tumor site
#' @return vcf_data_df
#'
#' @export
append_cancer_association_ranks <-
  function(var_data_df = NULL,
           ref_data = NULL,
           primary_site = 'Any') {

    if (any(grepl(paste0("^ENTREZGENE$"), names(var_data_df))) &
        any(grepl(paste0("^ENSEMBL_GENE_ID$"), names(var_data_df)))) {

      tissue_gene_ranks <- ref_data[['gene']][['otp_rank']] |>
        dplyr::select(
          c("ENTREZGENE", "PRIMARY_SITE", "TISSUE_ASSOC_RANK")) |>
        dplyr::filter(.data$PRIMARY_SITE == primary_site) |>
        dplyr::distinct()

      global_gene_ranks <- ref_data[['gene']][['otp_rank']] |>
        dplyr::select(c("ENTREZGENE", "GLOBAL_ASSOC_RANK")) |>
        dplyr::distinct()


      var_data_df_1 <- var_data_df |>
        dplyr::filter(!is.na(.data$ENTREZGENE))
      var_data_df_2 <- var_data_df |>
        dplyr::filter(
          is.na(.data$ENTREZGENE) &
            !is.na(.data$ENSEMBL_GENE_ID))
      var_data_df_3 <- var_data_df |>
        dplyr::filter(
          is.na(.data$ENTREZGENE) &
            is.na(.data$ENSEMBL_GENE_ID))

      if (NROW(var_data_df_1) > 0) {

        ## Add gene ranks (Open Targets Platform)
        ## - according to primary tumor types/sites
        ## - globally (across all tumor types/sites)
        var_data_df_1 <- var_data_df_1 |>
          dplyr::left_join(
            global_gene_ranks, by = "ENTREZGENE") |>
          dplyr::mutate(GLOBAL_ASSOC_RANK = dplyr::if_else(
            is.na(.data$GLOBAL_ASSOC_RANK),
            as.numeric(0),
            as.numeric(.data$GLOBAL_ASSOC_RANK)
          ))
        if (NROW(tissue_gene_ranks) > 0) {
          tissue_gene_ranks$PRIMARY_SITE <- NULL
          var_data_df_1 <- var_data_df_1 |>
            dplyr::left_join(
              tissue_gene_ranks, by = "ENTREZGENE") |>
            dplyr::mutate(TISSUE_ASSOC_RANK = dplyr::if_else(
              is.na(.data$TISSUE_ASSOC_RANK),
              as.numeric(0),
              as.numeric(.data$TISSUE_ASSOC_RANK)
            ))
        }else{
          var_data_df_1 <- var_data_df_1 |>
            dplyr::mutate(TISSUE_ASSOC_RANK = as.numeric(0))
        }

      }

      var_data_df <- dplyr::bind_rows(
        var_data_df_1,
        var_data_df_2,
        var_data_df_3) |>
        dplyr::mutate(GLOBAL_ASSOC_RANK = dplyr::if_else(
          is.na(.data$GLOBAL_ASSOC_RANK),
          as.numeric(0),
          as.numeric(.data$GLOBAL_ASSOC_RANK)
        )) |>
        dplyr::mutate(TISSUE_ASSOC_RANK = dplyr::if_else(
          is.na(.data$TISSUE_ASSOC_RANK),
          as.numeric(0),
          as.numeric(.data$TISSUE_ASSOC_RANK)
        )) |>
        dplyr::distinct()

    }

    return(var_data_df)

  }


#' Function that appends cancer gene evidence links
#'
#' @param var_data_df Data frame of sample variants from VCF
#' @param ref_data PCGR/CPSR reference data object
#' @return var_data_df
#'
#' @export
append_cancer_gene_evidence <-
  function(var_data_df = NULL,
           ref_data = NULL) {

    if (any(grepl(paste0("^ENTREZGENE$"), names(var_data_df))) &
        any(grepl(paste0("^ENSEMBL_GENE_ID$"), names(var_data_df)))) {

      var_data_df_1 <- var_data_df |>
        dplyr::filter(!is.na(.data$ENTREZGENE))
      var_data_df_2 <- var_data_df |>
        dplyr::filter(
          is.na(.data$ENTREZGENE) &
            !is.na(.data$ENSEMBL_GENE_ID))
      var_data_df_3 <- var_data_df |>
        dplyr::filter(
          is.na(.data$ENTREZGENE) &
            is.na(.data$ENSEMBL_GENE_ID))

      if (NROW(var_data_df_1) > 0) {

        var_data_df_1 <- var_data_df_1 |>
          dplyr::left_join(
            dplyr::filter(
              dplyr::select(
                ref_data[["gene"]][["gene_xref"]],
                c("ENTREZGENE",
                  "ENSEMBL_GENE_ID",
                  "CANCERGENE_EVIDENCE")),
              !is.na(.data$ENTREZGENE)),
            by = c("ENTREZGENE" = "ENTREZGENE",
                   "ENSEMBL_GENE_ID" = "ENSEMBL_GENE_ID")) |>
          dplyr::distinct()
      }

      if (NROW(var_data_df_2) > 0) {

        var_data_df_2 <- var_data_df_2 |>
          dplyr::left_join(
            dplyr::filter(
              dplyr::select(
                ref_data[["gene"]][["gene_xref"]],
                c("ENSEMBL_GENE_ID","CANCERGENE_EVIDENCE")),
              !is.na(.data$ENSEMBL_GENE_ID)),
            by = c("ENSEMBL_GENE_ID")) |>
          dplyr::distinct()

      }

      var_data_df <- dplyr::bind_rows(
        var_data_df_1,
        var_data_df_2,
        var_data_df_3) |>
        dplyr::distinct() |>
        dplyr::mutate(CANCERGENE_EVIDENCE = dplyr::if_else(
          .data$CANCERGENE_EVIDENCE == ".",
          as.character(NA),
          as.character(.data$CANCERGENE_EVIDENCE)
        ))
    }

    if("CGC_TIER" %in% colnames(var_data_df) &
       "CGC_GERMLINE" %in% colnames(var_data_df) &
       "CGC_SOMATIC" %in% colnames(var_data_df)){

      var_data_df$CANCER_GENE_CENSUS <- ""
      var_data_df <- var_data_df |>
        dplyr::mutate(
          CANCER_GENE_CENSUS = dplyr::case_when(
            !is.na(CGC_TIER) &
              CGC_GERMLINE == TRUE &
              CGC_SOMATIC == FALSE ~
              paste0("CGC tier ", CGC_TIER, " - Germline"),
            !is.na(CGC_TIER) &
              CGC_SOMATIC == TRUE &
              CGC_GERMLINE == FALSE ~
              paste0("CGC tier ", CGC_TIER, " - Somatic"),
            !is.na(CGC_TIER) &
              CGC_SOMATIC == TRUE &
              CGC_GERMLINE == TRUE ~
              paste0("CGC tier ", CGC_TIER, " - Somatic/Germline"),
            TRUE ~ as.character("")
          )
        )
      var_data_df$CGC_TIER <- NULL
      var_data_df$CGC_SOMATIC <- NULL
      var_data_df$CGC_GERMLINE <- NULL

    }

    return(var_data_df)

  }


#' A function that generates a HTML link for selected
#' identifiers (DBSNP, COSMIC, CLINVAR, ENTREZ)
#'
#' @param var_df data frame
#' @param vardb type of database
#' @param vartype 'snv_indel' or 'cna'
#' @param group_by_var variable used for grouping (VAR_ID)
#' @param url_prefix url prefix for link generation
#' @param link_key_var variable used in url for linking
#' @param link_display_var variable used in url for display
#' @return df_annotation_links
#'
#' @export
generate_annotation_link <- function(
    var_df,
    vardb = "DBSNP",
    vartype = "snv_indel",
    group_by_var = c("VAR_ID","ENTREZGENE"),
    url_prefix = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=",
    link_key_var = "DBSNP_RSID",
    link_display_var = "DBSNP_RSID") {

  df_annotation_links <- data.frame()
  invisible(assertthat::assert_that(
    is.data.frame(var_df),
    msg = paste0("Object 'var_df' must be of type data.frame"))
  )
  assertable::assert_colnames(
    var_df, c(group_by_var,
              link_key_var,
              link_display_var),
    only_colnames = F, quiet = T)

  if (length(
    unique(
      group_by_var %in% colnames(var_df)
    )) == 1 &
    link_key_var %in% colnames(var_df) &
    link_display_var %in% colnames(var_df)) {
    selected_vars <- unique(c(group_by_var,
                              link_key_var,
                              link_display_var))
    tmp_df <- var_df[, selected_vars]
    tmp_df <- tmp_df[!is.na(tmp_df[, link_key_var]) &
                       !is.na(tmp_df[, link_display_var]), ] |>
      dplyr::distinct()


    if (nrow(tmp_df) > 0) {
      df_annotation_links <- data.frame()
      #if (vardb == "DBSNP") {
      df_annotation_links <- as.data.frame(
        tmp_df |>
          tidyr::separate_rows(
            !!rlang::sym(link_key_var),
            sep = "&|,") |>
          dplyr::mutate(
            tmp_link = paste0(
              "<a href='", url_prefix,
              !!rlang::sym(link_key_var), "' target='_blank'>",
              !!rlang::sym(link_display_var), "</a>")
          )
      )

      if (nrow(df_annotation_links) > 0) {
        if(vartype == "cna"){
          df_annotation_links <- as.data.frame(
            df_annotation_links |>
              dplyr::group_by(
                .data$VAR_ID,
                .data$ENTREZGENE,
                .data$TRANSCRIPT_START,
                .data$TRANSCRIPT_END) |>
              dplyr::summarise(
                link = unlist(
                  paste(
                    .data$tmp_link, collapse = ", ")),
                .groups = "drop")
          )
        }else{
          if(vartype == "exp"){
            df_annotation_links <- as.data.frame(
              df_annotation_links |>
                dplyr::group_by(
                  .data$ENTREZGENE) |>
                dplyr::summarise(
                  link = unlist(
                    paste(
                      .data$tmp_link, collapse = ", ")),
                  .groups = "drop")
            )
          }
          if(vartype == "snv_indel"){
            df_annotation_links <- as.data.frame(
              df_annotation_links |>
                dplyr::group_by(
                  .data$VAR_ID,
                  .data$ENTREZGENE) |>
                dplyr::summarise(
                  link = unlist(
                    paste(
                      .data$tmp_link, collapse = ", ")),
                  .groups = "drop")
            )
          }
        }
      }

    }
  }
  else{
    cat("WARNING: Could not generate HTML URL links", sep = "\n")
    cat("- missing url_key or grouping varable in df", sep = "\n")
  }
  return(df_annotation_links)
}

