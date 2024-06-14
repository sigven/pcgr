#' Function that generates mutational signatures data for PCGR report
#'
#' @param variant_set Somatic callset (SNV)
#' @param vstats Variant statistics object
#' @param ref_data PCGR reference data object
#' @param settings PCGR run/configuration settings
#' @param n_bootstrap_iterations Number of bootstrap iterations for signature analysis
#' @param sig_contribution_cutoff Minimal signature contribution for reporting
#'
#' @export
generate_report_data_signatures <-
  function(variant_set = NULL,
           vstats = NULL,
           ref_data = NULL,
           settings = NULL,
           n_bootstrap_iterations = 200,
           sig_contribution_cutoff = 0.01) {

    n_snvs_required <- 30

    pcg_report_signatures <-
      pcgrr::init_m_signature_content()

    if(!is.null(variant_set) &
       !is.null(vstats) &
       !is.null(ref_data) &
       !is.null(settings)){
      pcgrr::log4r_info("------")
      pcgrr::log4r_info("Identifying mutational signatures")
    }else{
      pcgrr::log4r_warn("Missing input data for mutational signature analysis")
      return(pcg_report_signatures)
    }

    if("n_snv" %in% names(vstats)){
      if(vstats$n_snv < n_snvs_required){
        pcgrr::log4r_warn(
          paste0("Too few SNVs detected in sample (n = ",
                 vstats$n_snv,")",
                 " - omitting mutational signature analysis"))
        pcg_report_signatures$missing_data <- TRUE
        return(pcg_report_signatures)
      }
    }

    sig_settings <- settings$conf$somatic_snv$mutational_signatures

    cosmic_metadata <-
      ref_data$metadata |>
      dplyr::filter(.data$source_abbreviation == "cosmic_mutsigs") |>
      dplyr::select(c("source_version")) |>
      dplyr::mutate(
        source_version = stringr::str_replace_all(
          .data$source_version, "[\r\n]" , ""))

    pcgrr::log4r_info("------")
    pcgrr::log4r_info(
      paste0("Identifying weighted contributions of reference ",
             "mutational signatures (COSMIC ",
             cosmic_metadata$source_version,")"))
    assay <- tolower(settings$conf$assay_properties$type)

    vcf_name_mutsig_analysis <-
      file.path(settings$output_dir,
                paste(
                  settings$sample_id,
                  stringi::stri_rand_strings(
                    1, 15, pattern = "[A-Za-z0-9]"),
                  "mutational_patterns_input.all.vcf",
                  sep="."))

    pcgrr::write_processed_vcf(
      calls = variant_set,
      sample_name = settings$sample_id,
      output_directory = settings$output_dir,
      vcf_fname = vcf_name_mutsig_analysis,
      snv_only = F)

    fit_signatures_to_ttype <- !as.logical(
      sig_settings$all_reference_signatures
    )

    ## Retrieve relevant signatures for the tumor in question


    sites_with_sig_prevalence <-
      unique(ref_data$misc$mutational_signature[["PRIMARY_SITE"]])

    site_has_prevalence_data <- T
    if(!(settings$conf$sample_properties$site %in%
         sites_with_sig_prevalence)){
      site_has_prevalence_data <- F
      pcgrr::log4r_warn(
        paste0("No signature prevalence data available for site '",
               settings$conf$sample_properties$site,
               "' - considering all signatures for analysis"))
    }

    prevalent_site_signatures <- NULL
    if (fit_signatures_to_ttype == T & site_has_prevalence_data == T) {
      prevalent_site_signatures <-
        pcgrr::get_prevalent_site_signatures(
          site = settings$conf$sample_properties$site,
          min_prevalence_pct =
            as.numeric(sig_settings$prevalence_reference_signatures),
          ref_data = ref_data,
          incl_poss_artifacts =
            sig_settings$include_artefact_signatures)
    }else{
      prevalent_site_signatures <-
        pcgrr::get_prevalent_site_signatures(
          site = "Any",
          min_prevalence_pct =
            as.numeric(sig_settings$prevalence_reference_signatures),
          ref_data = ref_data,
          incl_poss_artifacts =
            sig_settings$include_artefact_signatures)
    }

    ## read MutationalPattern VCF file
    if (file.exists(glue::glue("{vcf_name_mutsig_analysis}.gz"))) {
      grl <- suppressMessages(suppressWarnings(
        MutationalPatterns::read_vcfs_as_granges(
          vcf_files = glue::glue("{vcf_name_mutsig_analysis}.gz"),
          sample_names = settings$sample_id,
          type = "all",
          genome = ref_data$assembly$bsg,
          predefined_dbs_mbs = T)
      )
      )
      if(class(grl)[1] == "CompressedGRangesList"){

        snv_grl <- suppressMessages(
          MutationalPatterns::get_mut_type(
            grl, type = "snv", predefined_dbs_mbs = T))
        indel_grl <- suppressMessages(
          MutationalPatterns::get_mut_type(
            grl, type = "indel", predefined_dbs_mbs = T))
        indel_counts <- NULL
        if(length(indel_grl[[1]]) > 0){
          indel_grl <- MutationalPatterns::get_indel_context(
            indel_grl, ref_genome = ref_data$assembly$bsg)
          indel_counts <- MutationalPatterns::count_indel_contexts(
            indel_grl)
        }

        type_occurrences <- suppressWarnings(
          MutationalPatterns::mut_type_occurrences(
            snv_grl, ref_genome = ref_data$assembly$bsg)
        )

        if(is.data.frame(type_occurrences)){
          if("T>A" %in% colnames(type_occurrences)){
            type_occurrences <- type_occurrences |>
              dplyr::rename("A>T" = "T>A")
          }
          if("T>C" %in% colnames(type_occurrences)){
            type_occurrences <- type_occurrences |>
              dplyr::rename("A>G" = "T>C")
          }
          if("T>G" %in% colnames(type_occurrences)){
            type_occurrences <- type_occurrences |>
              dplyr::rename("A>C" = "T>G")
          }
        }

        num_snvs_sig_analysis <- as.character(
          formattable::comma(length(snv_grl[[1]]), digits = 0))

        pcgrr::log4r_info(paste0("Number of SNVs for signature analysis: ",
                                 num_snvs_sig_analysis))

        pcg_report_signatures[["result"]][["indel_counts"]] <-
          indel_counts
        pcg_report_signatures[["result"]][["type_occurrences"]] <-
          type_occurrences
        pcg_report_signatures[["eval"]] <- TRUE


        ## assign variants to variant set
        pcg_report_signatures[["variant_set"]][["all"]] <-
          data.frame('VAR_ID' = rownames(S4Vectors::mcols(snv_grl[[1]])),
                     stringsAsFactors = F) |>
          tidyr::separate(.data$VAR_ID, c('CHROM', 'pos_ref_alt'),
                          sep=":", remove = T) |>
          tidyr::separate(.data$pos_ref_alt, c("POS","ref_alt"),
                          sep="_", remove = T) |>
          tidyr::separate(.data$ref_alt, c("REF","ALT"),
                          sep = "/", remove = T) |>
          dplyr::mutate(POS = as.integer(.data$POS))

        ## get context matrix
        mut_mat <-
          MutationalPatterns::mut_matrix(
            vcf_list = snv_grl,
            ref_genome = ref_data$assembly$bsg,
            extension = 1)
        pcg_report_signatures[["result"]][["mut_mat"]] <- mut_mat
        mut_mat <- mut_mat + 0.0001

        # get reference signatures (COSMIC v3.4)
        all_reference_signatures <-
          pcgrr::cosmic_sbs_signatures[['all']]

        if(as.logical(sig_settings$include_artefact_signatures) == FALSE){
          all_reference_signatures <-
            pcgrr::cosmic_sbs_signatures[['no_artefacts']]
        }

        if (length(snv_grl[[1]]) >= n_snvs_required) {

          sig_similarity <-
            MutationalPatterns::cos_sim_matrix(
              mut_mat, all_reference_signatures)

          if(length(sig_similarity) >= 67){
            pcg_report_signatures[["result"]][['signature_similarity']] <-
              tidyr::pivot_longer(
                as.data.frame(sig_similarity),
                names_to = "SIGNATURE_ID",
                values_to = "SIMILARITY",
                cols = dplyr::everything()) |>
              dplyr::arrange(
                dplyr::desc(.data$SIMILARITY)) |>
              dplyr::mutate(
                SIMILARITY = round(
                  .data$SIMILARITY, digits = 4)) |>
              dplyr::left_join(
                dplyr::select(
                  prevalent_site_signatures$aetiology,
                  c("SIGNATURE_ID", "AETIOLOGY_KEYWORD")),
                by = "SIGNATURE_ID"
              ) |>
              #dplyr::mutate(SITE_SPECIFIC = "NOT_DEFINED") |>
              dplyr::mutate(SITE_SPECIFIC = dplyr::if_else(
                as.logical(fit_signatures_to_ttype) == TRUE &
                  site_has_prevalence_data == TRUE,
                "NO",
                as.character("NOT_DEFINED")
              )) |>
              dplyr::mutate(SITE_SPECIFIC = dplyr::case_when(
                SIGNATURE_ID %in%
                  unique(prevalent_site_signatures$aetiology$SIGNATURE_ID) &
                  as.logical(fit_signatures_to_ttype) == TRUE &
                  as.logical(site_has_prevalence_data) == TRUE ~ "YES",
                TRUE ~ as.character(SITE_SPECIFIC)))
          }
        }

        if (length(snv_grl[[1]]) >= sig_settings[["mutation_limit"]]) {
          ## select subset of signatures based on those prevalent in tumor type/tissue
          selected_sigs <- intersect(
            colnames(all_reference_signatures),
            unique(prevalent_site_signatures$aetiology$SIGNATURE_ID)
          )
          selected_reference_signatures <-
            all_reference_signatures[, selected_sigs]

          ## reconstruct mutation profile from reference mutational signatures
          ## using bootstrapping
          n_bootstrap_iterations <- 200
          i <- 1
          bootstrap_data <- list()
          bootstrap_data[['goodness_of_fit']] <- data.frame()
          bootstrap_data[['contributions']] <- data.frame()
          while(i <= n_bootstrap_iterations){
            mut_mat_bs <- mutpat_resample_mut_mat(
              mut_matrix = mut_mat)
            fit_ref_bs <-
              MutationalPatterns::fit_to_signatures(
                mut_mat_bs, selected_reference_signatures)
            if(!is.null(fit_ref_bs)){
              b <- as.data.frame(stats::setNames(
                reshape2::melt(fit_ref_bs[["contribution"]]),
                c("SIGNATURE_ID", "sample_id",
                  "contribution_raw"))) |>
                dplyr::mutate(
                  prop = .data$contribution_raw / sum(.data$contribution_raw)) |>
                dplyr::mutate(bootstrap_iteration = i)
              bootstrap_data[['contributions']] <- dplyr::bind_rows(
                bootstrap_data[['contributions']], b)

              sim_original_reconstructed <-
                as.data.frame(
                  MutationalPatterns::cos_sim_matrix(
                    mut_mat_bs, fit_ref_bs[["reconstructed"]]))
              colnames(sim_original_reconstructed) <- c("cosine_sim")
              rownames(sim_original_reconstructed) <- NULL
              #magrittr::set_colnames("cosine_sim") |>
              #magrittr::set_rownames(NULL) |>
              sim_original_reconstructed <-
                sim_original_reconstructed |>
                dplyr::mutate(bootstrap_iteration = i,
                              sample_id = settings$sample_id)
              bootstrap_data[['goodness_of_fit']] <-
                dplyr::bind_rows(
                  bootstrap_data[['goodness_of_fit']],
                  sim_original_reconstructed)
              #i <- n_bootstrap_iterations + 1
            }
            i <- i + 1
          }

          ## calculate confidence intervals for each signature
          sig_prop_data <- data.frame()
          for(sig in unique(bootstrap_data[['contributions']]$SIGNATURE_ID)){
            sdata <- bootstrap_data[['contributions']] |>
              dplyr::filter(.data$SIGNATURE_ID == sig)
            ci_data <- stats::t.test(sdata$prop, conf.level = 0.95)
            prop_ci_lower <- ci_data$conf.int[1]
            prop_ci_upper <- ci_data$conf.int[2]

            sig_prop_data <- dplyr::bind_rows(
              sig_prop_data,
              data.frame(SIGNATURE_ID = sig,
                         prop_signature = mean(sdata$prop),
                         prop_signature_ci_lower = prop_ci_lower,
                         prop_signature_ci_upper = prop_ci_upper)
            )

          }

          ci_data_gof <- stats::t.test(
            bootstrap_data[['goodness_of_fit']]$cosine_sim,
            conf.level = 0.95)
          gof <- list()
          gof[['ci_lower']] <- ci_data_gof$conf.int[1]
          gof[['ci_upper']] <- ci_data_gof$conf.int[2]
          gof[['estimate']] <- mean(bootstrap_data[['goodness_of_fit']]$cosine_sim)

          contributions_per_signature <-
            sig_prop_data |>
            dplyr::mutate(SIGNATURE_ID = as.character(.data$SIGNATURE_ID)) |>
            dplyr::mutate(sample_id = settings$sample_id,
                          n_bs_iterations = n_bootstrap_iterations) |>
            dplyr::select(c("SIGNATURE_ID",
                            "sample_id",
                            "n_bs_iterations",
                            "prop_signature",
                            "prop_signature_ci_lower",
                            "prop_signature_ci_upper")) |>
            dplyr::filter(.data$prop_signature > sig_contribution_cutoff) |>
            dplyr::arrange(dplyr::desc(.data$prop_signature)) |>
            dplyr::left_join(
              dplyr::select(
                ref_data$misc$mutational_signature,
                c("SIGNATURE_ID",
                  "AETIOLOGY",
                  "COMMENTS",
                  "AETIOLOGY_KEYWORD")),
              by = c("SIGNATURE_ID")) |>
            dplyr::rename(
              group = "AETIOLOGY_KEYWORD",
              signature_id = "SIGNATURE_ID",
              aetiology = "AETIOLOGY",
              comments = "COMMENTS") |>
            dplyr::mutate(
              contribution =
                paste0(
                  round(.data$prop_signature * 100, digits = 2), "%")) |>
            dplyr::distinct()

          contributions_per_group <- as.data.frame(
            contributions_per_signature |>
              dplyr::group_by(.data$group) |>
              dplyr::summarise(
                prop_group = sum(.data$prop_signature),
                signature_id_group = paste(
                  sort(.data$signature_id), collapse=", "),
                .groups = "drop") |>
              dplyr::arrange(
                dplyr::desc(.data$prop_group)) |>
              dplyr::mutate(group = dplyr::if_else(
                .data$prop_group > 0.05,
                .data$group,
                "Other")) |>
              dplyr::group_by(.data$group) |>
              dplyr::summarise(
                prop_group = sum(.data$prop_group),
                signature_id_group = paste(
                  sort(.data$signature_id_group), collapse=", "),
                .groups = "drop") |>
              dplyr::mutate(prop_group = round(
                .data$prop_group, digits = 3)) |>
              dplyr::arrange(
                dplyr::desc(.data$prop_group))

          )

          cols <- contributions_per_signature |>
            dplyr::arrange(dplyr::desc(.data$prop_signature)) |>
            dplyr::select(c("signature_id")) |>
            dplyr::distinct() |>
            utils::head(25)

          color_vec <- utils::head(
            pcgrr::color_palette[["tier"]][["values"]],
            min(25, nrow(cols)))

          names(color_vec) <- cols$signature_id
          color_vec2 <- color_vec
          names(color_vec2) <- NULL

          cols <- cols |>
            dplyr::mutate(col = color_vec2)
          contributions_per_signature <- contributions_per_signature |>
            dplyr::left_join(cols, by = "signature_id")

          ## emit warning if more than 25 different signatures are found
          ## choose only signatures attributed to 25 different signatures
          missing_signatures <- contributions_per_signature |>
            dplyr::filter(is.na(.data$col))
          if (nrow(missing_signatures) > 0) {
            log4r_warn(paste0("Found contributions from more than 25 signatures - ",
                              "showing signatures from 25 different signatures only"))
            contributions_per_signature <- contributions_per_signature |>
              dplyr::filter(!is.na(.data$col))

            contributions_per_group <- as.data.frame(
              contributions_per_group |>
                tidyr::separate_rows(
                  .data$signature_id_group, sep = ", ") |>
                dplyr::anti_join(missing_signatures,
                                 by = c("signature_id_group" = "signature_id")) |>
                dplyr::group_by(
                  c("group","prop_group")
                ) |>
                dplyr::summarise(
                  signature_id_group = paste(
                    sort(.data$signature_id_group), collapse=", "),
                  .groups = "drop"
                )
            )
          }


          contributions <- list()
          contributions[["per_group"]] <-  contributions_per_group
          contributions[["per_signature"]] <-  contributions_per_signature
          tsv_data <- data.frame()

          ## Get output for tab-separated file
          ## - contribution per signature id and reference signatures used
          if (!is.null(prevalent_site_signatures$aetiology) &
              NROW(contributions[["per_signature"]]) > 0) {
            if ("SIGNATURE_ID" %in% colnames(prevalent_site_signatures$aetiology)) {
              reference_sigs <- paste(sort(prevalent_site_signatures$aetiology$SIGNATURE_ID),
                                      collapse=",")
              tsv_data <- contributions[["per_signature"]] |>
                pcgrr::remove_cols_from_df(
                  cnames = c("contribution","col","AETIOLOGY","COMMENTS")) |>
                dplyr::mutate(
                  all_reference_signatures = !fit_signatures_to_ttype,
                  tumor_type = settings$conf$sample_properties$site,
                  reference_collection = "COSMIC_v34",
                  reference_signatures = reference_sigs,
                  fitting_accuracy =
                    round(gof$estimate * 100, digits = 1)) |>
                dplyr::select(c("sample_id"), dplyr::everything())
            }
          }

          vr <- grl[[settings$sample_id]]
          GenomeInfoDb::seqlengths(vr) <-
            GenomeInfoDb::seqlengths(
              ref_data$assembly$bsg)[GenomeInfoDb::seqlevels(ref_data$assembly$bsg) %in% unique(GenomeInfoDb::seqlevels(vr))]
          chromosomes <- utils::head(
            GenomeInfoDb::seqnames(ref_data$assembly$bsg), 24)

          pcg_report_signatures[["result"]][["vr"]] <- vr
          pcg_report_signatures[["result"]][["chromosomes"]] <- chromosomes
          pcg_report_signatures[["result"]][["contributions"]] <- contributions
          pcg_report_signatures[["result"]][["tsv"]] <- tsv_data
          pcg_report_signatures[["result"]][["no_site_prevalence"]] <-
            !site_has_prevalence_data
          pcg_report_signatures[["result"]][["reference_data"]] <-
            prevalent_site_signatures$aetiology
          pcg_report_signatures[["result"]][["scale_fill_values"]] <- color_vec
          pcg_report_signatures[["result"]][["scale_fill_names"]] <-
            names(color_vec)
          pcg_report_signatures[["result"]][["goodness_of_fit"]] <-
            gof
        }else{
          pcg_report_signatures[["missing_data"]] <- TRUE
          if (length(snv_grl[[1]]) > 0) {
            pcgrr::log4r_info(
              paste0("Too few SNVs (n = ",
                     nrow(pcg_report_signatures[["variant_set"]][["all"]]),
                     ") for reconstruction of mutational signatures by ",
                     "MutationalPatterns, limit set to ",
                     sig_settings$mutation_limit))
          }
        }
      }
    }

    fs::file_delete(glue::glue("{vcf_name_mutsig_analysis}.gz"))

    return(pcg_report_signatures)
  }


#' Function that retrieves prevalent signatures for a given tumor type/primary site
#' Data is collected from COSMIC v3.4.
#'
#' @param site Primary tumor site
#' @param custom_collection Custom collection of signatures from COSMIC
#' @param ref_data PCGR reference data object
#' @param min_prevalence_pct Minimum prevalence (pct) of signature in
#' cohorts associated with primary site -
#' used to select reference signatures for inclusion in signature reconstruction
#' @param incl_poss_artifacts logical indicating if artefact signatures
#' are to be included
#'
#' @export
get_prevalent_site_signatures <-
  function(site = "Any",
           custom_collection = NULL,
           ref_data = NULL,
           min_prevalence_pct = 0.1,
           incl_poss_artifacts = T) {

    cosmic_metadata <-
      ref_data$metadata |>
      dplyr::filter(
        .data$source_abbreviation == "cosmic_mutsigs") |>
      dplyr::select(c("source_version")) |>
      dplyr::mutate(
        source_version = stringr::str_replace_all(
          .data$source_version, "[\r\n]" , ""))

    if (is.null(custom_collection)) {
      pcgrr::log4r_info(paste0(
        "Retrieving prevalent (prevalence >= ",
        min_prevalence_pct, " percent) reference signatures for ",
        site, ", using COSMIC ",
        cosmic_metadata$source_version,
        " collection"))
    }
    pcgrr::log4r_info(paste0(
      "Inclusion of mutational signature artefacts (e.g. sequencing artefacts): ",
      as.logical(incl_poss_artifacts)))

    invisible(
      assertthat::assert_that(
        !is.null(ref_data$misc$mutational_signature),
        msg =
          paste0(
            "Cannot load ref. aetiologies (COSMIC ",
            cosmic_metadata$source_version,
            ") of mutational signatures")))
    invisible(
      assertthat::assert_that(
        is.data.frame(ref_data$misc$mutational_signature),
        msg = "Reference aetiologies must be of type data.frame()"))
    invisible(
      assertthat::assert_that(
        min_prevalence_pct >= 0.1 &
          min_prevalence_pct <= 20,
        msg = "Argument 'min_prevalence_pct' must be more than 0.1 and less than 20"))

    valid_signature_ids <-
      unique(ref_data$misc$mutational_signature$SIGNATURE_ID)
    signatures_prevalence <- data.frame()

    if (!is.null(custom_collection)) {
      invisible(
        assertthat::assert_that(
          is.character(custom_collection),
          msg = "Argument 'custom_collection' must be a character vector"))

      pcgrr::log4r_info(paste0(
        "Retrieving reference signatures from COSMIC ",
        cosmic_metadata$source_version,
        " collection based on user-defined collection (",
        paste(unique(custom_collection), collapse=", "), ")")
      )
      i <- 1
      while(i <= length(custom_collection)) {
        if (!(custom_collection[i] %in% valid_signature_ids)) {
          log4r_warn(paste0(
            "Could not find specified custom signature id  '",
            custom_collection[i], "' in COSMIC ",
            cosmic_metadata$source_version,
            " reference collection",
            " - ignoring"))
        }
        i <- i + 1
      }

      signatures_prevalence <-
        ref_data$misc$mutational_signature |>
        dplyr::select(c("SIGNATURE_ID",
                      "AETIOLOGY_KEYWORD",
                      "AETIOLOGY",
                      "ASSOCIATED_SIGNATURES",
                      "COMMENTS")) |>
        dplyr::filter(.data$SIGNATURE_ID %in% custom_collection) |>
        dplyr::distinct()

    }else{

      unique_sites_with_signature_prevalence <-
        unique(ref_data$misc$mutational_signature[["PRIMARY_SITE"]])

      ## No primary site defined - 'Any'
      if (!(site %in% unique_sites_with_signature_prevalence)) {
        signatures_prevalence <-
          ref_data$misc$mutational_signature |>
          dplyr::select(c("SIGNATURE_ID",
                          "AETIOLOGY_KEYWORD",
                          "AETIOLOGY",
                          "ASSOCIATED_SIGNATURES",
                          "COMMENTS")) |>
          dplyr::distinct()
      }else{
        signatures_prevalence <-
          ref_data$misc$mutational_signature |>
          dplyr::filter(.data$PRIMARY_SITE == site) |>
          dplyr::select(.data$SIGNATURE_ID,
                        .data$PRIMARY_SITE,
                        .data$PREVALENCE_PCT,
                        .data$AETIOLOGY_KEYWORD,
                        .data$AETIOLOGY,
                        .data$ASSOCIATED_SIGNATURES,
                        .data$COMMENTS) |>
          dplyr::distinct()

        if (min_prevalence_pct > 0.1) {
          signatures_prevalence <- signatures_prevalence |>
            dplyr::filter(!is.na(.data$PREVALENCE_PCT)) |>
            dplyr::filter(.data$PREVALENCE_PCT >= min_prevalence_pct)
        }
        signatures_prevalence <- signatures_prevalence |>
          dplyr::distinct() |>
          dplyr::arrange(dplyr::desc(.data$PREVALENCE_PCT)) |>
          dplyr::select(-c("PREVALENCE_PCT","PRIMARY_SITE"))
      }
    }

    if (incl_poss_artifacts == F) {
      signatures_prevalence <- signatures_prevalence |>
        dplyr::filter(!stringr::str_detect(
          .data$AETIOLOGY_KEYWORD,"artefact"))
    }
    signatures_prevalence <- signatures_prevalence |>
      dplyr::distinct()

    ## Subset signature matrix - keeping only columns (signatures)
    ## to those defined by primary site/custom collection
    sigs <- unique(signatures_prevalence$SIGNATURE_ID)
    #pcgrr::log4r_info(paste0("Limiting reference collection to signatures: ",
    #                          paste(sigs, collapse = ", ")))

    result <- list("aetiology" = signatures_prevalence)

    return(result)
  }

#' Function that generates data for rainfall plot (mutation density
#' along genome, considering SNVs only)
#'
#' @param variant_set data frame with SNVs/InDels (must contain "CHROM",
#' "POS","REF","ALT")
#' @param colors character vector of six color codes (denoting the
#' different mutation types)
#' @param autosomes logical indicating if plotting should only
#' consider autosomes
#' @param build genome assembly (grch37/grch38)
#'
#' @export
generate_report_data_rainfall <- function(variant_set,
                                          colors = NULL,
                                          autosomes = FALSE,
                                          build = NULL) {

  pcg_report_rainfall <- pcgrr::init_rainfall_content()
  if (NROW(variant_set) == 0) {
    return(pcg_report_rainfall)
  }

  invisible(assertthat::assert_that
            (assertthat::is.flag(autosomes),
              msg = "Argument 'autosomes' should be TRUE/FALSE"))
  invisible(assertthat::assert_that(
    is.data.frame(variant_set),
    msg = paste0("Argument variant_set needs be of type data.frame")))
  assertable::assert_colnames(
    variant_set, c("CHROM", "REF", "ALT",
                   "POS", "VARIANT_CLASS"), only_colnames = F, quiet = T)
  invisible(
    assertthat::assert_that(
      build == "grch37" | build == "grch38",
      msg = paste0("Value for argument build ('", build,
                   "') not allowed, available reference build values are:",
                   "'grch37' or 'grch38'")))

  pcgrr::log4r_info("------")
  pcgrr::log4r_info(paste0("Calculating data for rainfall plot"))


  sbs_types <- c("C>T", "A>G", "A>C", "A>T", "C>G", "C>A")
  if (is.null(colors)) {
    colors <- utils::head(pcgrr::color_palette$tier$values, 6)
  }else{
    invisible(
      assertthat::assert_that(
        length(colors) == 6,
        msg = "'colors' should be a character vector of six color codes"))
  }

  chr_prefix <- FALSE
  chromosome_names <- unique(variant_set[, "CHROM"])
  for (m in chromosome_names) {
    if (startsWith(m, "chr")) {
      chr_prefix <- TRUE
    }
  }
  if (chr_prefix == F) {
    variant_set <- variant_set |>
      dplyr::mutate(CHROM = paste0("chr", .data$CHROM))
  }

  dat <- dplyr::select(
    variant_set,
    c("CHROM",
      "POS",
      "REF",
      "ALT",
      "VARIANT_CLASS")) |>
    dplyr::filter(.data$VARIANT_CLASS == "SNV")

  if (nrow(dat) < 10 | nrow(dat) > 30000 | length(chromosome_names) < 2) {
    pcgrr::log4r_info(
      paste0("Too few variants (< 10) and chromosomes ",
             " represented (< 2) OR too many variants ",
             "( > 30,000) - skipping rainfall plot"))
    pcg_report_rainfall[["eval"]] <- F
  }else{
    dat <- dat |>
      pcgrr::assign_mutation_type() |>
      dplyr::mutate(
        MUTATION_TYPE =
          dplyr::if_else(stringr::str_detect(.data$MUTATION_TYPE, "^(C|A)>"),
                         stringr::str_replace(.data$MUTATION_TYPE,
                                              ":[A-Z]>[A-Z]$", ""),
                         as.character(.data$MUTATION_TYPE))) |>
      # dplyr::mutate(
      #   MUTATION_TYPE =
      #     dplyr::if_else(stringr::str_detect(.data$MUTATION_TYPE, "^A>"),
      #                    stringr::str_replace(.data$MUTATION_TYPE,
      #                                         "^[A-Z]>[A-Z]:", ""),
      #                    as.character(.data$MUTATION_TYPE))) |>
      pcgrr::sort_chromosomal_segments()

    bsg <- get_genome_obj(build)
    chr_length <- utils::head(GenomeInfoDb::seqlengths(bsg), 24)
    chromosomes <- utils::head(GenomeInfoDb::seqnames(bsg), 24)
    if (autosomes == TRUE) {
      chr_length <- utils::head(chr_length, 22)
    }

    # cumulative sum of chromosome lengths
    chr_cum <- c(0, cumsum(as.numeric(chr_length)))

    # Plot chromosome labels without "chr"
    names(chr_cum) <- names(chr_length)
    labels <- gsub("chr", "", names(chr_length))

    # position of chromosome labels
    m <- c()
    for (i in 2:length(chr_cum))
      m <- c(m, (chr_cum[i - 1] + chr_cum[i]) / 2)

    # mutation characteristics
    type <- c()
    loc <- c()
    dist <- c()
    chrom <- c()

    # for each chromosome
    #chromosomes <-
    for (i in 1:length(chromosomes)) {
      chr_subset <- dplyr::filter(dat, .data$CHROM == chromosomes[i])
      n <- nrow(chr_subset)
      if (n <= 1) {
        next
      }
      type <- c(type, chr_subset$MUTATION_TYPE[-1])
      loc <- c(loc, (chr_subset$POS + chr_cum[i])[-1])
      dist <- c(dist, diff(chr_subset$POS))
      chrom <- c(chrom, rep(chromosomes[i], n - 1))
    }

    invisible(assertthat::assert_that(length(type) == length(loc) &
                            length(loc) == length(dist) &
                            length(chrom) == length(dist),
                            msg = "Length of type/loc/dist not identical"))
    data <- data.frame(type = type,
                      location = loc,
                      distance = dist,
                      chromosome = chrom,
                      stringsAsFactors = F)

    # Removes colors based on missing mutation types.  This prevents colors from
    # shifting when comparing samples with low mutation counts.
    typesin <- sbs_types %in% sort(unique(data$type))
    colors_selected <- colors[typesin]
    ylim <- 1e+09
    pcg_report_rainfall[["eval"]] <- T
    pcg_report_rainfall[["rfdata"]] <-
      list("data" = data, "intercept" = m, "ylim" = ylim,
           "chr_cum" = chr_cum, "colors" = colors_selected,
           "labels" = labels, "cex" = 0.8, "cex_text" = 3)
  }
  return(pcg_report_rainfall)
}

#' MutationalPatterns Resample a mutation matrix
#'
#' Since this is an unexported function, copied verbatim from
#' https://github.com/UMCUGenetics/MutationalPatterns/blob/ca9caf/R/fit_to_signatures_bootstrapped.R#L161
#'
#' The mutation matrix is resampled per column (sample).
#' Resampling is done with replacement using the row weights as propabilities.
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#'
#' @return A resamples mutation matrix
#'
#' @noRd
#'
mutpat_resample_mut_mat <- function(mut_matrix) {
  mut_mat_resampled <- apply(mut_matrix, 2, function(x) {
    total_muts <- sum(x)
    sample_weights <- x / total_muts
    feature_rows <- sample(seq_along(x), total_muts, replace = TRUE, prob = sample_weights)
    row_counts <- table(feature_rows)
    index <- as.numeric(names(row_counts))
    x[index] <- as.vector(row_counts)
    x[-index] <- 0
    return(x)
  })
  return(mut_mat_resampled)
}
