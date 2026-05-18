#' Function that computes various variant statistics from a data frame
#' with variant records
#'
#' @param var_df data frame with variants
#' @param pct_other_limit numeric value specifying the percentage limit
#' for the 'Other' category
#'
#' @export
#'
stats_type_snv_indel <- function(
    var_df = NULL, pct_other_limit = 4) {

  assertthat::assert_that(
    !is.null(var_df),
    is.data.frame(var_df),
    msg = "Argument 'var_df' must be a valid data.frame"
  )

  assertable::assert_colnames(
    var_df, c("VARIANT_CLASS", "CONSEQUENCE","CODING_STATUS"),
    only_colnames = F, quiet = T
  )

  consequence_stats <-
    var_df |>
    dplyr::mutate(CONSEQUENCE = stringr::str_replace_all(
      .data$CONSEQUENCE, "(, [0-9A-Za-z_]{1,}) {1,}$",""
    )) |>
    dplyr::group_by(.data$CONSEQUENCE) |>
    dplyr::summarise(
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$N))

  if (NROW(consequence_stats) > 5) {
    consequence_stats_top <- utils::head(consequence_stats, 4)
    consequence_stats_other <- consequence_stats |>
      dplyr::slice_tail(n = -4) |>
      dplyr::summarise(
        N = sum(.data$N),
        CONSEQUENCE = "other_consequences"
      )
    consequence_stats <- dplyr::bind_rows(
      consequence_stats_top, consequence_stats_other) |>
      dplyr::arrange(dplyr::desc(.data$N))
  }

  consequence_stats <- consequence_stats |>
    dplyr::mutate(Pct = .data$N / sum(.data$N) * 100)

  consequence_stats_coding <-
    var_df |>
    dplyr::filter(.data$CODING_STATUS == "coding")

  if (NROW(consequence_stats_coding) > 0) {
    consequence_stats_coding <-
      consequence_stats_coding |>
      dplyr::mutate(CONSEQUENCE = stringr::str_replace_all(
        .data$CONSEQUENCE, "(, [0-9A-Za-z_]{1,}) {1,}$",""
      )) |>
      dplyr::group_by(.data$CONSEQUENCE) |>
      dplyr::summarise(
        N = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(.data$N))

    if (NROW(consequence_stats_coding) > 5) {
      consequence_stats_coding_top <- utils::head(consequence_stats_coding, 4)
      consequence_stats_coding_other <- consequence_stats_coding |>
        dplyr::slice_tail(n = -4) |>
        dplyr::summarise(
          N = sum(.data$N),
          CONSEQUENCE = "other_consequences"
        )
      consequence_stats_coding <- dplyr::bind_rows(
        consequence_stats_coding_top,
        consequence_stats_coding_other) |>
        dplyr::arrange(dplyr::desc(.data$N))
    }

    consequence_stats_coding <-
      consequence_stats_coding |>
      dplyr::mutate(Pct = .data$N / sum(.data$N) * 100)
  }


  variant_class_stats <-
    var_df |>
    dplyr::group_by(.data$VARIANT_CLASS) |>
    dplyr::summarise(
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(Pct = .data$N / sum(.data$N) * 100) |>
    dplyr::arrange(dplyr::desc(.data$Pct))

  coding_stats <-
    var_df |>
    dplyr::group_by(.data$CODING_STATUS) |>
    dplyr::summarise(
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(Pct = .data$N / sum(.data$N) * 100) |>
    dplyr::arrange(dplyr::desc(.data$Pct))

  result <- list()
  result[['consequence']] <- consequence_stats
  result[['consequence_coding']] <- consequence_stats_coding
  result[['variant_class']] <- variant_class_stats
  result[['coding']] <- coding_stats

  return(result)
}

#' Function that computes various variant statistics for
#' SNVs/InDels from a callset object
#'
#' @param callset list object with callset data (SNVs/InDels)
#'
#' @export
#'
stats_report_snv_indel <- function(
    callset = NULL) {

  invisible(assertthat::assert_that(
    !is.null(callset),
    is.list(callset),
    msg = "Argument 'callset' can not be NULL and must be a list object"
  ))

  ## check that variant is a member of callset and is a
  ## data frame with the required columns for variant statistics

  invisible(assertthat::assert_that(
    "variant" %in% names(callset) &
      is.data.frame(callset$variant),
    msg = "Argument 'callset' must contain a data frame member named 'variant'"
  ))

  assertable::assert_colnames(
    callset$variant,
    c("VARIANT_CLASS",
      "CODING_STATUS",
      "ACTIONABILITY_TIER",
      "TUMOR_SUPPRESSOR",
      "ONCOGENE",
      "gnomADe_AF",
      "TARGETED_INHIBITORS_ALL2"),
    only_colnames = F, quiet = T
  )

  vstats <- pcgrr::init_vstats_snv_indel()

  vstats[["n_var"]] <-
    callset$variant |>
    NROW()
  if (vstats[["n_var"]] > 0) {
    vstats[["n_var_snv"]] <-
      callset$variant |>
      dplyr::filter(
        .data$VARIANT_CLASS == "SNV") |>
      NROW()
    vstats[["n_var_indel"]] <-
      callset$variant |>
      dplyr::filter(
        .data$VARIANT_CLASS == "insertion" |
          .data$VARIANT_CLASS == "deletion" |
          .data$VARIANT_CLASS == "indel") |>
      nrow()
    vstats[["n_var_sub"]] <-
      callset$variant |>
      dplyr::filter(
        .data$VARIANT_CLASS == "substituion") |>
      nrow()
    vstats[["n_var_coding"]] <-
      callset$variant |>
      dplyr::filter(.data$CODING_STATUS == "coding") |>
      nrow()
    vstats[["n_var_noncoding"]] <-
      callset$variant |>
      dplyr::filter(.data$CODING_STATUS == "noncoding") |>
      nrow()
    vstats[["n_var_exonic"]] <-
      callset$variant |>
      dplyr::filter(.data$EXONIC_STATUS == "exonic") |>
      NROW()
    vstats[["n_var_nonexonic"]] <-
      callset$variant |>
      dplyr::filter(.data$EXONIC_STATUS == "nonexonic") |>
      NROW()
  }

  if ("bm_evidence" %in% names(callset)) {
    for (class in names(pcgrr::bm_categories)) {
      if (class %in% names(callset$bm_evidence)) {

        if("classification" %in% names(callset$bm_evidence[[class]]) &
           "eitems" %in% names(callset$bm_evidence[[class]])) {
          var_classification <-
            callset$bm_evidence[[class]][['classification']]
          eitems <-
            callset$bm_evidence[[class]][['eitems']]

          if(NROW(var_classification) == 0 |
             NROW(eitems) == 0){
            pcgrr::log4r_debug(
              paste0("No classification or evidence items for class '",
                     class,"' - skipping statistics for this class"))
            next
          }

          ## print debug message
          pcgrr::log4r_debug(
            paste0("Processing biomarker evidence for class '",
                   class,"' with ",
                   NROW(var_classification),
                   " classified variants and ",
                   NROW(eitems),
                   " evidence items")
          )

          invisible(assertthat::assert_that(
            is.data.frame(var_classification) &
              is.data.frame(eitems),
            msg = paste0("Classification and eitems for class '",
                         class,"' must be data frames")
          ))

          invisible(assertable::assert_colnames(
            var_classification,
            c("VAR_ID",
              "VARIANT_CLASS",
              "ACTIONABILITY_TIER",
              "ENTREZGENE"),
            only_colnames = F, quiet = T
          ))

          invisible(assertable::assert_colnames(
            eitems,
            c("VAR_ID",
              "BM_ACTIONABILITY_SUPPORT",
              "ENTREZGENE"),
            only_colnames = F, quiet = T
          ))

          if (NROW(var_classification) > 0) {
            var_classification <- var_classification |>
              dplyr::select(c("VAR_ID",
                              "VARIANT_CLASS",
                              "ACTIONABILITY_TIER",
                              "ENTREZGENE")) |>
              dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
              dplyr::distinct()

            for (tier in 1:6) {

              c1 <- class
              c2 <- paste0("tier",tier)

              tier_vars <-
                var_classification |>
                dplyr::filter(.data$ACTIONABILITY_TIER == tier) |>
                dplyr::distinct()

              if (NROW(tier_vars) > 0) {
                num_tier_genes <- tier_vars |>
                  dplyr::select(
                    "ENTREZGENE") |>
                  dplyr::distinct() |>
                  NROW()

                vstats[[c1]][[c2]][['n_eitems_tier_defining']] <-
                  eitems |>
                  dplyr::filter(
                    .data$BM_ACTIONABILITY_SUPPORT == "tier-defining" &
                      .data$VAR_ID %in% tier_vars$VAR_ID &
                      .data$ENTREZGENE %in% tier_vars$ENTREZGENE) |>
                  NROW()

                vstats[[c1]][[c2]][['n_eitems_additional']] <-
                  eitems |>
                  dplyr::filter(
                    .data$BM_ACTIONABILITY_SUPPORT == "additional" &
                      .data$VAR_ID %in% tier_vars$VAR_ID &
                      .data$ENTREZGENE %in% tier_vars$ENTREZGENE) |>
                  NROW()

                vstats[[c1]][[c2]][['n_gene']] <-
                  num_tier_genes
                vstats[[c1]][[c2]][['n_var']] <-
                  NROW(tier_vars)

                if (tier == 3){

                  vstats[[c1]][[c2]][['n_var_tsg']] <-
                    tier_vars |>
                    dplyr::select("VAR_ID",
                                  "ENTREZGENE",
                                  "ACTIONABILITY_TIER",
                                  "VARIANT_CLASS") |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "ACTIONABILITY_TIER",
                             "VARIANT_CLASS",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                        !is.na(.data$TUMOR_SUPPRESSOR) &
                        .data$TUMOR_SUPPRESSOR == TRUE &
                        .data$ONCOGENE == FALSE) |>
                    NROW()

                  vstats[[c1]][[c2]][['n_var_oncg']] <-
                    tier_vars |>
                    dplyr::select("VAR_ID",
                                  "ENTREZGENE",
                                  "ACTIONABILITY_TIER",
                                  "VARIANT_CLASS") |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "ACTIONABILITY_TIER",
                             "VARIANT_CLASS",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                        !is.na(.data$ONCOGENE) &
                        .data$ONCOGENE == TRUE &
                        .data$TUMOR_SUPPRESSOR == FALSE) |>
                    NROW()

                  vstats[[c1]][[c2]][['n_var_dual_role']] <-
                    tier_vars |>
                    dplyr::select(
                      "VAR_ID",
                      "ENTREZGENE",
                      "ACTIONABILITY_TIER",
                      "VARIANT_CLASS") |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "VARIANT_CLASS",
                             "ACTIONABILITY_TIER",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                       .data$TUMOR_SUPPRESSOR == TRUE &
                        .data$ONCOGENE == TRUE) |>
                    NROW()

                }
              }
            }
          }
        }
      }
    }
  }

  return(vstats)
}


#' Function that computes various variant statistics for CNAs
#' from a callset object
#'
#' @param callset list object with callset data (CNA)
#'
#' @export
stats_report_cna <- function(
    callset = NULL) {

  invisible(assertthat::assert_that(
    !is.null(callset),
    is.list(callset),
    msg = "Argument 'callset' can not be NULL and must be a list object"
  ))

  ## check that variant is a member of callset and is a
  ## data frame with the required columns for variant statistics

  invisible(assertthat::assert_that(
    "variant" %in% names(callset) &
      is.data.frame(callset$variant),
    msg = "Argument 'callset' must contain a data frame member named 'variant'"
  ))

  assertable::assert_colnames(
    callset$variant,
    c("VARIANT_CLASS",
      "VAR_ID",
      "ENTREZGENE",
      "LOH",
      "ACTIONABILITY_TIER",
      "TUMOR_SUPPRESSOR",
      "ONCOGENE",
      "TARGETED_INHIBITORS_ALL2"),
    only_colnames = F, quiet = T
  )

  vstats <- pcgrr::init_vstats_cna()

  for (cat in c("tsg","oncogene","drugtarget")) {

    df <- callset$variant |>
      dplyr::select(
        c("VAR_ID",
          "ENTREZGENE",
          "VARIANT_CLASS",
          "LOH",
          "ACTIONABILITY_TIER",
          "TUMOR_SUPPRESSOR",
          "ONCOGENE",
          "TARGETED_INHIBITORS_ALL2")) |>
      dplyr::filter(
        !is.na(.data$VARIANT_CLASS) &
          .data$VARIANT_CLASS %in%
          c("homdel",
            "hemdel",
            "hetdel",
            "neutral",
            "gain",
            "amplification")) |>
      dplyr::distinct()

    if (cat == "tsg") {
      df <- df |>
        dplyr::filter(
          !is.na(.data$TUMOR_SUPPRESSOR) &
            .data$TUMOR_SUPPRESSOR == TRUE)
    }
    if (cat == "oncogene") {
      df <- df |>
        dplyr::filter(
          !is.na(.data$ONCOGENE) &
            .data$ONCOGENE == TRUE)
    }
    if (cat == "drugtarget") {
      df <- df |>
        dplyr::filter(
          !is.na(
            .data$TARGETED_INHIBITORS_ALL2))
    }

    for (stat in c("n_gene_homdel",
                   "n_gene_hetdel",
                   "n_gene_gain",
                   "n_gene_amplification",
                   "n_gene_LOH_deletion",
                   "n_gene_LOH_amplification",
                   "n_gene_LOH_copy_neutral")) {

      num_hits <- 0
      if (stat == "n_gene_homdel") {
        num_hits <- df |>
          dplyr::filter(
            (.data$VARIANT_CLASS == "homdel" |
               .data$VARIANT_CLASS == "hemdel")) |>
          NROW()
      }
      if (stat == "n_gene_hetdel") {
        num_hits <- df |>
          dplyr::filter(
            .data$VARIANT_CLASS == "hetdel") |>
          NROW()
      }
      if (stat == "n_gene_amplification") {
        num_hits <- df |>
          dplyr::filter(
            .data$VARIANT_CLASS == "amplification") |>
          NROW()
      }
      if (stat == "n_gene_gain") {
        num_hits <- df |>
          dplyr::filter(
            .data$VARIANT_CLASS == "gain") |>
          NROW()
      }
      if (stat == "n_gene_LOH_deletion") {
        num_hits <- df |>
          dplyr::filter(
            !is.na(.data$LOH) &
              .data$LOH == "deletion") |>
          NROW()
      }
      if (stat == "n_gene_LOH_amplification") {
        num_hits <- df |>
          dplyr::filter(
            !is.na(.data$LOH) &
              .data$LOH == "amplification") |>
          NROW()
      }
      if (stat == "n_gene_LOH_copy_neutral") {
        num_hits <- df |>
          dplyr::filter(
            !is.na(.data$LOH) &
              .data$LOH == "copy_neutral") |>
          NROW()
      }

      vstats[[cat]][[stat]] <- num_hits
    }
  }

  if ("bm_evidence" %in% names(callset)) {
    for (class in names(pcgrr::bm_categories)) {
      if (class %in% names(callset$bm_evidence)) {

        if("classification" %in% names(callset$bm_evidence[[class]]) &
           "eitems" %in% names(callset$bm_evidence[[class]])) {
          var_classification <-
            callset$bm_evidence[[class]][['classification']]
          eitems <-
            callset$bm_evidence[[class]][['eitems']]

          if(NROW(var_classification) == 0 |
             NROW(eitems) == 0){
            pcgrr::log4r_debug(
              paste0("No classification or evidence items for class '",
                     class,"' - skipping statistics for this class"))
            next
          }

          invisible(assertthat::assert_that(
            is.data.frame(var_classification) &
              is.data.frame(eitems),
            msg = paste0("Classification and eitems for class '",
                         class,"' must be data frames")
          ))

          invisible(assertable::assert_colnames(
            var_classification,
            c("VAR_ID",
              "VARIANT_CLASS",
              "ACTIONABILITY_TIER",
              "ENTREZGENE"),
            only_colnames = F, quiet = T
          ))

          invisible(assertable::assert_colnames(
            eitems,
            c("VAR_ID",
              "BM_ACTIONABILITY_SUPPORT",
              "ENTREZGENE"),
            only_colnames = F, quiet = T
          ))

          if (NROW(var_classification) > 0) {
            var_classification <- var_classification |>
              dplyr::select(c("VAR_ID",
                              "VARIANT_CLASS",
                              "ACTIONABILITY_TIER",
                              "ENTREZGENE")) |>
              dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
              dplyr::distinct()

            for (tier in 1:6) {

              c1 <- class
              c2 <- paste0("tier",tier)

              tier_vars <-
                var_classification |>
                dplyr::filter(.data$ACTIONABILITY_TIER == tier) |>
                dplyr::distinct()

              if (NROW(tier_vars) > 0) {
                num_tier_genes <- tier_vars |>
                  dplyr::select(
                    "ENTREZGENE","VAR_ID") |>
                  dplyr::distinct() |>
                  NROW()

                vstats[[c1]][[c2]][['n_eitems_tier_defining']] <-
                  eitems |>
                  dplyr::filter(
                    .data$BM_ACTIONABILITY_SUPPORT == "tier-defining" &
                      .data$VAR_ID %in% tier_vars$VAR_ID &
                      .data$ENTREZGENE %in% tier_vars$ENTREZGENE) |>
                  NROW()

                vstats[[c1]][[c2]][['n_eitems_additional']] <-
                  eitems |>
                  dplyr::filter(
                    .data$BM_ACTIONABILITY_SUPPORT == "additional" &
                      .data$VAR_ID %in% tier_vars$VAR_ID &
                      .data$ENTREZGENE %in% tier_vars$ENTREZGENE) |>
                  NROW()

                vstats[[c1]][[c2]][['n_gene']] <-
                  num_tier_genes

                if (tier == 3){

                  vstats[[c1]][[c2]][['n_gene_tsg_homdel']] <-
                    tier_vars |>
                    dplyr::select("VAR_ID",
                                  "ENTREZGENE",
                                  "ACTIONABILITY_TIER",
                                  "VARIANT_CLASS") |>
                    dplyr::distinct() |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "ACTIONABILITY_TIER",
                             "VARIANT_CLASS",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                      (.data$VARIANT_CLASS == "homdel" |
                         .data$VARIANT_CLASS == "hemdel") &
                        !is.na(.data$TUMOR_SUPPRESSOR) &
                        .data$TUMOR_SUPPRESSOR == TRUE) |>
                    NROW()

                  vstats[[c1]][[c2]][['n_gene_tsg_hetdel']] <-
                    tier_vars |>
                    dplyr::select("VAR_ID",
                                  "ENTREZGENE",
                                  "ACTIONABILITY_TIER",
                                  "VARIANT_CLASS") |>
                    dplyr::distinct() |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "ACTIONABILITY_TIER",
                             "VARIANT_CLASS",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                      .data$VARIANT_CLASS == "hetdel" &
                        !is.na(.data$TUMOR_SUPPRESSOR) &
                        .data$TUMOR_SUPPRESSOR == TRUE) |>
                    NROW()

                  vstats[[c1]][[c2]][['n_gene_oncg_amplification']] <-
                    tier_vars |>
                    dplyr::select("VAR_ID",
                                  "ENTREZGENE",
                                  "ACTIONABILITY_TIER",
                                  "VARIANT_CLASS") |>
                    dplyr::distinct() |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "VARIANT_CLASS",
                             "ACTIONABILITY_TIER",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                      .data$VARIANT_CLASS == "amplification" &
                        !is.na(.data$ONCOGENE) &
                        .data$ONCOGENE == TRUE) |>
                    NROW()

                  vstats[[c1]][[c2]][['n_gene_oncg_gain']] <-
                    tier_vars |>
                    dplyr::select("VAR_ID",
                                  "ENTREZGENE",
                                  "ACTIONABILITY_TIER",
                                  "VARIANT_CLASS") |>
                    dplyr::distinct() |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "VARIANT_CLASS",
                             "ACTIONABILITY_TIER",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                      .data$VARIANT_CLASS == "gain" &
                        !is.na(.data$ONCOGENE) &
                        .data$ONCOGENE == TRUE) |>
                    NROW()

                  vstats[[c1]][[c2]][['n_gene_drugtarget_amplification']] <-
                    tier_vars |>
                    dplyr::select(
                      "VAR_ID",
                      "ENTREZGENE",
                      "ACTIONABILITY_TIER",
                      "VARIANT_CLASS") |>
                    dplyr::distinct() |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "VARIANT_CLASS",
                             "ACTIONABILITY_TIER",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                      .data$VARIANT_CLASS == "amplification" &
                        .data$ONCOGENE == FALSE &
                        !is.na(.data$TARGETED_INHIBITORS_ALL2)) |>
                    NROW()

                }
              }
            }
          }
        }
      }
    }
  }
  return(vstats)
}

#' Function that computes various variant statistics for fusions
#' from a callset object
#'
#' @param callset list object with callset data (fusions)
#'
#' @export
#'
stats_report_fusion <- function(
    callset = NULL) {

  vstats <- pcgrr::init_vstats_fusion()

  invisible(assertthat::assert_that(
    !is.null(callset),
    is.list(callset),
    msg = "Argument 'callset' can not be NULL and must be a list object"
  ))

  ## check that variant is a member of callset and is a
  ## data frame with the required columns for variant statistics

  invisible(assertthat::assert_that(
    "variant" %in% names(callset) &
      is.data.frame(callset$variant),
    msg = "Argument 'callset' must contain a data frame member named 'variant'"
  ))

  assertable::assert_colnames(
    callset$variant,
    c("VARIANT_CLASS",
      "ACTIONABILITY_TIER",
      "ENTREZGENE",
      "ONCOGENE_5P",
      "ONCOGENE_3P",
      "TARGETED_INHIBITORS_ALL2_5P",
      "TARGETED_INHIBITORS_ALL2_3P"),
    only_colnames = F, quiet = T
  )

  vstats[["n_fusion"]] <-
    callset$variant |>
    NROW()
  if (vstats[["n_fusion"]] > 0) {
    vstats[["n_fusion_oncg"]] <-
      callset$variant |>
      dplyr::filter(
        .data$ONCOGENE_5P == TRUE |
          .data$ONCOGENE_3P == TRUE) |>
      NROW()
    vstats[["n_fusion_oncg_5P"]] <-
      callset$variant |>
      dplyr::filter(
        .data$ONCOGENE_5P == TRUE |
          .data$ONCOGENE_3P == TRUE) |>
      NROW()
    vstats[["n_fusion_oncg_3P"]] <-
      callset$variant |>
      dplyr::filter(
        .data$ONCOGENE_5P == TRUE |
          .data$ONCOGENE_3P == TRUE) |>
      NROW()

  }

  if ("bm_evidence" %in% names(callset)) {
    for (class in names(pcgrr::bm_categories)) {
      if (class %in% names(callset$bm_evidence)) {

        if("classification" %in% names(callset$bm_evidence[[class]]) &
           "eitems" %in% names(callset$bm_evidence[[class]])) {
          var_classification <-
            callset$bm_evidence[[class]][['classification']]
          eitems <-
            callset$bm_evidence[[class]][['eitems']]

          if(NROW(var_classification) == 0 |
             NROW(eitems) == 0) {
            next
          }

          ## print debug message
          pcgrr::log4r_debug(
            paste0("Processing biomarker evidence for class '",
                   class,"' with ",
                   NROW(var_classification),
                   " classified variants and ",
                   NROW(eitems),
                   " evidence items")
          )

          invisible(assertthat::assert_that(
            is.data.frame(var_classification) &
              is.data.frame(eitems),
            msg = paste0("Classification and eitems for class '",
                         class,"' must be data frames")
          ))

          invisible(assertable::assert_colnames(
            var_classification,
            c("VAR_ID",
              "VARIANT_CLASS",
              "ACTIONABILITY_TIER",
              "ENTREZGENE"),
            only_colnames = F, quiet = T
          ))

          invisible(assertable::assert_colnames(
            eitems,
            c("VAR_ID",
              "BM_ACTIONABILITY_SUPPORT",
              "ENTREZGENE"),
            only_colnames = F, quiet = T
          ))

          if (NROW(var_classification) > 0) {
            var_classification <- var_classification |>
              dplyr::select(c("VAR_ID",
                              "VARIANT_CLASS",
                              "ACTIONABILITY_TIER",
                              "ENTREZGENE")) |>
              dplyr::filter(!is.na(.data$ACTIONABILITY_TIER)) |>
              dplyr::distinct()

            for (tier in 1:6) {

              c1 <- class
              c2 <- paste0("tier",tier)

              tier_vars <-
                var_classification |>
                dplyr::filter(.data$ACTIONABILITY_TIER == tier) |>
                dplyr::distinct()

              if (NROW(tier_vars) > 0) {
                num_tier_fusions <- tier_vars |>
                  dplyr::select(
                    "ENTREZGENE") |>
                  dplyr::distinct() |>
                  NROW()

                vstats[[c1]][[c2]][['n_eitems_tier_defining']] <-
                  eitems |>
                  dplyr::filter(
                    .data$BM_ACTIONABILITY_SUPPORT == "tier-defining" &
                      .data$VAR_ID %in% tier_vars$VAR_ID &
                      .data$ENTREZGENE %in% tier_vars$ENTREZGENE) |>
                  NROW()

                vstats[[c1]][[c2]][['n_eitems_additional']] <-
                  eitems |>
                  dplyr::filter(
                    .data$BM_ACTIONABILITY_SUPPORT == "additional" &
                      .data$VAR_ID %in% tier_vars$VAR_ID &
                      .data$ENTREZGENE %in% tier_vars$ENTREZGENE) |>
                  NROW()

                vstats[[c1]][[c2]][['n_fusion']] <-
                  num_tier_fusions

                if (tier == 3){

                  vstats[[c1]][[c2]][['n_fusion_oncg']] <-
                    tier_vars |>
                    dplyr::select("VAR_ID",
                                  "ENTREZGENE",
                                  "ACTIONABILITY_TIER",
                                  "VARIANT_CLASS") |>
                    dplyr::inner_join(
                      callset[["variant"]],
                      by = c("VAR_ID",
                             "ACTIONABILITY_TIER",
                             "VARIANT_CLASS",
                             "ENTREZGENE")
                    ) |>
                    dplyr::filter(
                      (!is.na(.data$ONCOGENE_3P) &
                        .data$ONCOGENE_3P == TRUE) |
                        (!is.na(.data$ONCOGENE_5P) &
                           .data$ONCOGENE_5P == TRUE)) |>
                    NROW()

                }
              }
            }
          }
        }
      }
    }
  }

  return(vstats)
}

#' Function that generate stats for a germline variant callset,
#' including number of variants and number of variants with
#' evidence items for each BM evidence type
#'
#' @param var_df data.frame with germline variant data (SNVs/InDels)
#' @param vartype type of variant ('snv_indel')
#'
#' @export
stats_report_germline <- function(
    var_df = NULL,
    vartype = "snv_indel") {

  call_stats <- list()

  invisible(
    assertthat::assert_that(
      !is.null(var_df) &
        is.data.frame(var_df),
      msg = "Argument 'var_df' must be a data.frame()"
    )
  )

  if (vartype == 'snv_indel' &
     "CLASSIFICATION" %in% colnames(var_df)) {
    if ("BM_EVIDENCE_TYPE" %in% colnames(var_df)) {

      call_stats[['n_eitems_predictive']] <- var_df |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Predictive") |>
        NROW()
      call_stats[['n_eitems_prognostic']] <- var_df |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Prognostic") |>
        NROW()
      call_stats[['n_eitems_diagnostic']] <- var_df |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Diagnostic") |>
        NROW()
      call_stats[['n_eitems_predisposing']] <- var_df |>
        dplyr::filter(.data$BM_EVIDENCE_TYPE == "Predisposing") |>
        NROW()

      call_stats[['n_var_eitems']] <- var_df |>
        dplyr::filter(!is.na(.data$BM_EVIDENCE_TYPE)) |>
        dplyr::select("GENOMIC_CHANGE") |>
        dplyr::distinct() |>
        NROW()
    }

    if ( "VARIANT_CLASS" %in% colnames(var_df) &
         "CODING_STATUS" %in% colnames(var_df) &
         "GENOMIC_CHANGE" %in% colnames(var_df)){
      call_stats[['n']] <- var_df |>
        NROW()

      call_stats[['n_snv']] <- var_df |>
        dplyr::filter(.data$VARIANT_CLASS == "SNV") |>
        NROW()

      call_stats[['n_indel']] <- var_df |>
        dplyr::filter(.data$VARIANT_CLASS == "insertion" |
                        .data$VARIANT_CLASS == "deletion" |
                        .data$VARIANT_CLASS == "indel") |>
        NROW()

      call_stats[['n_coding']] <- var_df |>
        dplyr::filter(.data$CODING_STATUS == "coding") |>
        NROW()

      call_stats[['n_noncoding']] <- var_df |>
        dplyr::filter(.data$CODING_STATUS == "noncoding") |>
        NROW()

      call_stats[['n_p']] <- var_df |>
        dplyr::filter(.data$CLASSIFICATION == "Pathogenic") |>
        NROW()
      call_stats[['n_lp']] <- var_df |>
        dplyr::filter(.data$CLASSIFICATION == "Likely Pathogenic") |>
        NROW()
      call_stats[['n_vus']] <- var_df |>
        dplyr::filter(.data$CLASSIFICATION == "VUS") |>
        NROW()
      call_stats[['n_lb']] <- var_df |>
        dplyr::filter(.data$CLASSIFICATION == "Likely Benign") |>
        NROW()
      call_stats[['n_b']] <- var_df |>
        dplyr::filter(.data$CLASSIFICATION == "Benign") |>
        NROW()
    }

  }

  return(call_stats)
}

#' Function that initiates report element with actionable
#' variant statistics information
#'
#' @export
init_vstats_actionable <- function(vartype = "snv_indel") {
  vstats <- list()
  for (cat in names(pcgrr::bm_categories)) {
    vstats[[cat]] <- list()
    for (i in 1:6) {
      if(i > 3 & vartype != "snv_indel"){
        next
      }
      vstats[[cat]][[paste0('tier',i)]] <- list()
      for (stat in c('n_var',
                     'n_gene',
                     'n_fusion',
                     'n_eitems_tier_defining',
                     'n_eitems_additional')) {
        if(stat == "n_var" & vartype == "cna"){
          next
        }
        if((stat == "n_gene" | stat == "n_var") &
           vartype == "fusion"){
          next
        }
        if(stat == "n_fusion" & vartype != "fusion"){
          next
        }

        vstats[[cat]][[paste0('tier',i)]][[stat]] <- 0
        if (i == 3 & vartype == "snv_indel") {
          for (gene_role in c("n_var_tsg",
                              "n_var_oncg",
                              "n_var_dual_role")) {
            vstats[[cat]][[paste0('tier',i)]][[gene_role]] <- 0
          }
        }
        if (i == 3 & vartype == "fusion") {
          vstats[[cat]][[paste0('tier',i)]][["n_fusion_oncg"]] <- 0
        }
        if (i == 3 & vartype == "cna") {
          for (gene_role in c("n_gene_tsg_homdel",
                              "n_gene_tsg_hetdel",
                              "n_gene_oncg_gain",
                              "n_gene_oncg_amplification",
                              "n_gene_drugtarget_amplification")) {
            vstats[[cat]][[paste0('tier',i)]][[gene_role]] <- 0
          }
        }
      }
    }
  }
  return(vstats)
}


#' Function that initiates RNA fusion statistics
#'
#' @export
init_vstats_fusion <- function() {
  vstats <- init_vstats_actionable(vartype = "fusion")
   for (t in c("n_fusion",
               "n_fusion_oncg",
               "n_fusion_oncg_5P",
               "n_fusion_oncg_3P")) {
    vstats[[t]] <- 0
  }
  return(vstats)
}

#' Function that initiates CNA statistics
#'
#' @export
init_vstats_cna <- function() {

  vstats <- init_vstats_actionable(vartype = "cna")
  for (cat in c("tsg",
               "oncogene",
               "drugtarget",
               "segments")) {
    vstats[[cat]] <- list()
    stat_names <- c(
      "n_segment",
      "n_segment_homdel",
      "n_segment_hetdel",
      "n_segment_gain",
      "n_segment_amplification",
      "n_segment_LOH_deletion",
      "n_segment_LOH_amplification",
      "n_segment_LOH_copy_neutral"
    )
    if (cat != "segments") {
      stat_names <- c(
        "n_gene_homdel",
        "n_gene_hetdel",
        "n_gene_gain",
        "n_gene_amplification",
        "n_gene_LOH_deletion",
        "n_gene_LOH_amplification",
        "n_gene_LOH_copy_neutral"
      )
    }
    for (stat in stat_names) {
        vstats[[cat]][[stat]] <- 0
      }
  }
  return(vstats)
}

#' Function that initiates report element with SNV/InDel statistics information
#'
#' @export
init_vstats_snv_indel <- function() {

  vstats <- init_vstats_actionable(vartype = "snv_indel")
  for (t in c("n_var",
              "n_var_snv",
              "n_var_indel",
              "n_var_sub",
              "n_var_exonic",
              "n_var_nonexonic",
              "n_var_coding",
              "n_var_noncoding")) {
    vstats[[t]] <- 0
  }
  return(vstats)
}

