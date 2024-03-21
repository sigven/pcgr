#' Function that generates mutational signatures data for PCGR report
#'
#' @param variant_set Somatic callset (SNV)
#' @param ref_data PCGR reference data object
#' @param settings PCGR configuration settings object
#' @param n_bootstrap_iterations Number of bootstrap iterations for signature analysis
#' @param sig_contribution_cutoff Minimal signature contribution for reporting
#'
#' @export
generate_report_data_signatures <-
  function(variant_set = NULL,
           ref_data = NULL,
           settings = NULL,
           n_bootstrap_iterations = 200,
           sig_contribution_cutoff = 0) {


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
  #assay <- tolower(pcgr_config$assay_props$type)
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


  pcg_report_signatures <-
    pcgrr::init_m_signature_content()

  fit_signatures_to_ttype <- !as.logical(
    sig_settings$all_reference_signatures
  )

  ## Retrieve relevant signatures for the tumor in question
  prevalent_site_signatures <- NULL
  if (fit_signatures_to_ttype == T) {
    prevalent_site_signatures <-
      pcgrr::get_prevalent_site_signatures(
        site = settings$conf$sample_properties$site,
        min_prevalence_pct =
          sig_settings$prevalence_reference_signatures,
        ref_data = ref_data,
        incl_poss_artifacts =
          sig_settings$include_artefact_signatures)
  }else{
    prevalent_site_signatures <-
      pcgrr::get_prevalent_site_signatures(
        site = "Any",
        min_prevalence_pct =
          sig_settings$prevalence_reference_signatures,
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

      snv_grl <- suppressMessages(MutationalPatterns::get_mut_type(
        grl, type = "snv"))
      indel_grl <- suppressMessages(MutationalPatterns::get_mut_type(
        grl, type = "indel", predefined_dbs_mbs = T))
      indel_counts <- NULL
      if(length(indel_grl[[1]]) > 0){
        indel_grl <- MutationalPatterns::get_indel_context(
          indel_grl, ref_genome = ref_data$assembly$bsg)
        indel_counts <- MutationalPatterns::count_indel_contexts(indel_grl)
      }

      type_occurrences <- MutationalPatterns::mut_type_occurrences(
        snv_grl, ref_genome = ref_data$assembly$bsg)
      pcgrr::log4r_info(paste0("Number of SNVs for signature analysis: ",
                               length(snv_grl[[1]])))

      pcg_report_signatures[["eval"]] <- TRUE

      if (length(snv_grl[[1]]) >= sig_settings[["mutation_limit"]]) {

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
        mut_mat <- mut_mat + 0.0001

        ## get reference signatures (COSMIC v3.2)
        all_reference_signatures <-
          MutationalPatterns::get_known_signatures(
          muttype = "snv",
          genome = stringr::str_replace(
            ref_data$assembly$grch_name, "grc", "GRC"
          ),
          incl_poss_artifacts =
            as.logical(
              sig_settings$include_artefact_signatures
            )
        )

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
          mut_mat_bs <- MutationalPatterns:::.resample_mut_mat(
            mut_matrix = mut_mat)
          fit_ref_bs <-
            MutationalPatterns::fit_to_signatures(
              mut_mat_bs, selected_reference_signatures)
          if(!is.null(fit_ref_bs)){
            b <- as.data.frame(stats::setNames(
              reshape2::melt(fit_ref_bs[["contribution"]]),
              c("SIGNATURE_ID", "sample_id",
                "contribution_raw"))) |>
              dplyr::mutate(prop = contribution_raw / sum(contribution_raw)) |>
              dplyr::mutate(bootstrap_iteration = i)
            bootstrap_data[['contributions']] <- dplyr::bind_rows(
              bootstrap_data[['contributions']], b)

            sim_original_reconstructed <-
              as.data.frame(
                MutationalPatterns::cos_sim_matrix(
                  mut_mat_bs, fit_ref_bs[["reconstructed"]])) |>
              magrittr::set_colnames("cosine_sim") |>
              magrittr::set_rownames(NULL) |>
              dplyr::mutate(bootstrap_iteration = i,
                            sample_id = settings$sample_id)
            bootstrap_data[['goodness_of_fit']] <- dplyr::bind_rows(
              bootstrap_data[['goodness_of_fit']], sim_original_reconstructed)
            #i <- n_bootstrap_iterations + 1
          }
          i <- i + 1
        }

        ## calculate confidence intervals for each signature
        sig_prop_data <- data.frame()
        for(sig in unique(bootstrap_data[['contributions']]$SIGNATURE_ID)){
          sdata <- bootstrap_data[['contributions']] |>
            dplyr::filter(SIGNATURE_ID == sig)
          ci_data <- t.test(sdata$prop, conf.level = 0.95)
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

        ci_data_gof <- t.test(bootstrap_data[['goodness_of_fit']]$cosine_sim,
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
              paste0(round(.data$prop_signature * 100, digits = 2), "%")) |>
          dplyr::distinct()

        contributions_per_group <- as.data.frame(
          contributions_per_signature |>
            dplyr::group_by(.data$group) |>
            dplyr::summarise(
              prop_group = sum(.data$prop_signature),
              signature_id_group = paste(
                .data$signature_id, collapse=", "),
              .groups = "drop")

        )

        cols <- contributions_per_group |>
          dplyr::arrange(dplyr::desc(.data$prop_group)) |>
          dplyr::select(.data$group) |>
          dplyr::distinct() |>
          utils::head(25)

        color_vec <- utils::head(
          pcgrr::color_palette[["tier"]][["values"]], min(25, nrow(cols)))

        names(color_vec) <- cols$group
        color_vec2 <- color_vec
        names(color_vec2) <- NULL

        cols <- cols |>
          dplyr::mutate(col = color_vec2)
        contributions_per_signature <- contributions_per_signature |>
          dplyr::left_join(cols, by = "group")

        ## emit warning if more than 25 different aetiologies are found
        ## choose only signatures attributed to 25 different aetiologies
        missing_aetiologies <- contributions_per_signature |>
          dplyr::filter(is.na(.data$col))
        if (nrow(missing_aetiologies) > 0) {
          log4r_warn(paste0("Found contributions from more than 25 aetiologies - ",
                            "showing signatures from 25 different aetiologies only"))
          contributions_per_signature <- contributions_per_signature |>
            dplyr::filter(!is.na(.data$col))

          contributions_per_group <- contributions_per_group |>
            dplyr::anti_join(missing_aetiologies, by = "group")
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
                  round(gof$estimate * 100, digits = 1))
          }
        }

        vr <- grl[[settings$sample_id]]
        GenomeInfoDb::seqlengths(vr) <-
          GenomeInfoDb::seqlengths(ref_data$assembly$bsg)[GenomeInfoDb::seqlevels(ref_data$assembly$bsg) %in% unique(GenomeInfoDb::seqlevels(vr))]
        chromosomes <- utils::head(GenomeInfoDb::seqnames(ref_data$assembly$bsg), 24)

        pcg_report_signatures[["result"]][["vr"]] <- vr
        pcg_report_signatures[["result"]][["indel_counts"]] <- indel_counts
        pcg_report_signatures[["result"]][["type_occurrences"]] <- type_occurrences
        pcg_report_signatures[["result"]][["mut_mat"]] <- mut_mat
        pcg_report_signatures[["result"]][["chromosomes"]] <- chromosomes
        pcg_report_signatures[["result"]][["contributions"]] <- contributions
        pcg_report_signatures[["result"]][["tsv"]] <- tsv_data
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
          pcg_report_signatures[["variant_set"]][["all"]] <-
            as.data.frame(snv_grl[[1]]) |>
            dplyr::rename(POS = .data$start, CHROM = .data$seqnames) |>
            dplyr::select(.data$CHROM, .data$POS, .data$REF, .data$ALT) |>
            magrittr::set_rownames(NULL)
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

  system(glue::glue("rm -f {vcf_name_mutsig_analysis}*"))

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
           min_prevalence_pct = 5,
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
      incl_poss_artifacts))

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
        min_prevalence_pct == 1 |
          min_prevalence_pct == 2 |
          min_prevalence_pct == 5 |
          min_prevalence_pct == 10 |
          min_prevalence_pct == 15 |
          min_prevalence_pct == 20,
        msg = "Argument 'min_prevalence_pct' must be any of '0, 2, 5, 10, 15 or 20'"))

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
      if (!(site %in% unique_sites_with_signature_prevalence)) {
        pcgrr::log4r_info(
          paste0("Primary tumor site '", site, "' ",
                 "does not have any signatures with significant ",
                 "prevalence - considering all"))
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
                        .data$PREVALENCE_ABOVE_5PCT,
                        .data$PREVALENCE_ABOVE_10PCT,
                        .data$PREVALENCE_ABOVE_15PCT,
                        .data$PREVALENCE_ABOVE_20PCT,
                        .data$AETIOLOGY_KEYWORD,
                        .data$AETIOLOGY,
                        .data$ASSOCIATED_SIGNATURES,
                        .data$COMMENTS) |>
          dplyr::distinct()

        if (min_prevalence_pct > 0) {
          if (min_prevalence_pct == 5) {
            signatures_prevalence <- signatures_prevalence |>
              dplyr::filter(.data$PREVALENCE_ABOVE_5PCT == T |
                              is.na(.data$PREVALENCE_ABOVE_5PCT))
          }else if (min_prevalence_pct == 10) {
            signatures_prevalence <- signatures_prevalence |>
              dplyr::filter(.data$PREVALENCE_ABOVE_10PCT == T |
                              is.na(.data$PREVALENCE_ABOVE_10PCT))
          }
          else if (min_prevalence_pct == 15) {
            signatures_prevalence <- signatures_prevalence |>
              dplyr::filter(.data$PREVALENCE_ABOVE_15PCT == T |
                              is.na(.data$PREVALENCE_ABOVE_15PCT))
          }else if (min_prevalence_pct == 20) {
            signatures_prevalence <- signatures_prevalence |>
              dplyr::filter(.data$PREVALENCE_ABOVE_20PCT == T |
                              is.na(.data$PREVALENCE_ABOVE_20PCT))
          }else if (min_prevalence_pct == 2 | min_prevalence_pct == 1) {
            signatures_prevalence <- signatures_prevalence |>
              dplyr::filter(!is.na(.data$PREVALENCE_PCT)) |>
              dplyr::filter(.data$PREVALENCE_PCT >= min_prevalence_pct)
          }
        }
        signatures_prevalence <- signatures_prevalence |>
          dplyr::select(-c(.data$PRIMARY_SITE,
                           .data$PREVALENCE_ABOVE_5PCT,
                           .data$PREVALENCE_ABOVE_10PCT,
                           .data$PREVALENCE_ABOVE_15PCT,
                           .data$PREVALENCE_ABOVE_20PCT)) |>
          dplyr::distinct() |>
          dplyr::arrange(dplyr::desc(.data$PREVALENCE_PCT)) |>
          dplyr::select(-.data$PREVALENCE_PCT)
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


  sbs_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (is.null(colors)) {
    colors <- head(pcgrr::color_palette$tier$values, 6)
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

  dat <- dplyr::select(variant_set, .data$CHROM, .data$POS,
                       .data$REF, .data$ALT, .data$VARIANT_CLASS) |>
    dplyr::filter(.data$VARIANT_CLASS == "SNV")

  if (nrow(dat) < 10 | nrow(dat) > 10000 | length(chromosome_names) < 2) {
    pcgrr::log4r_info(
      paste0("Too few variants (< 10) and chromosomes ",
             " represented (< 2) OR too many variants ",
             "( > 10,000) - skipping rainfall plot"))
    pcg_report_rainfall[["eval"]] <- F
  }else{
    dat <- dat |>
      pcgrr::assign_mutation_type() |>
      dplyr::mutate(
        MUTATION_TYPE =
          dplyr::if_else(stringr::str_detect(.data$MUTATION_TYPE, "^C>"),
                         stringr::str_replace(.data$MUTATION_TYPE,
                                              ":[A-Z]>[A-Z]$", ""),
                         as.character(.data$MUTATION_TYPE))) |>
      dplyr::mutate(
        MUTATION_TYPE =
          dplyr::if_else(stringr::str_detect(.data$MUTATION_TYPE, "^A>"),
                         stringr::str_replace(.data$MUTATION_TYPE,
                                              "^[A-Z]>[A-Z]:", ""),
                         as.character(.data$MUTATION_TYPE))) |>
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
