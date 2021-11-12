#' Function that generates mutational signatures data for PCGR report
#'
#' @param vcf_fname VCF file processed with PCGR annotation pipeline -
#' possibly filtered for depth/allelic fraction
#' @param pcgr_data object with PCGR annotation data
#' @param sample_name sample identifier
#' @param pcgr_config Object with PCGR configuration parameters
#' @param type_specific logical indicating if all reference signatures are to be
#' included (F) rather than those known to be prevalent in the tumor (T)
#'
#' @export
generate_report_data_signatures_mp <-
  function(vcf_fname,
           pcgr_data,
           sample_name,
           pcgr_config,
           type_specific = T) {

  log4r_info("------")
  log4r_info(paste0("Identifying weighted contributions of reference ",
                    "mutational signatures (COSMIC v3.2) using ",
                    "MutationalPatterns"))
  assay <- tolower(pcgr_config$assay_props$type)


  pcg_report_signatures <-
    pcgrr::init_report(config = pcgr_config,
                       class = "m_signature_mp")

  ## Retrieve relevant signatures for the tumor in question
  prevalent_site_signatures <- NULL
  if(type_specific == T){
    prevalent_site_signatures <-
      pcgrr::get_prevalent_site_signatures(
        site = pcgr_config[["t_props"]][["tumor_type"]],
        pcgr_data = pcgr_data,
        incl_poss_artifacts =
          pcgr_config[["msigs"]][["include_artefact_signatures"]])
  }
  if(type_specific == F){
    prevalent_site_signatures <-
      pcgrr::get_prevalent_site_signatures(
        site = "Any",
        pcgr_data = pcgr_data,
        incl_poss_artifacts =
          pcgr_config[["msigs"]][["include_artefact_signatures"]])
  }

  ## read MutationalPattern VCF file
  if(file.exists(vcf_fname)){
    vcfs <- suppressWarnings(
      MutationalPatterns::read_vcfs_as_granges(
        vcf_files = vcf_fname,
        sample_names = sample_name,
        genome = pcgr_data[["assembly"]][["ref_genome"]],
        predefined_dbs_mbs = T),
      )

    log4r_info(paste0("Number of SNVs for signature analysis: ",
                             length(vcfs[[1]])))

    pcg_report_signatures[["eval"]] <- TRUE

    if (length(vcfs[[1]]) >= pcgr_config[["msigs"]][["mutation_limit"]]) {

      ## assign variants to variant set
      pcg_report_signatures[["variant_set"]][["all"]] <-
        data.frame('VAR_ID' = rownames(mcols(vcfs[[1]])),
                        stringsAsFactors = F) %>%
        tidyr::separate(.data$VAR_ID, c('CHROM', 'pos_ref_alt'),
                        sep=":", remove = T) %>%
        tidyr::separate(.data$pos_ref_alt, c("POS","ref_alt"),
                        sep="_", remove = T) %>%
        tidyr::separate(.data$ref_alt, c("REF","ALT"),
                        sep = "/", remove = T) %>%
        dplyr::mutate(POS = as.integer(.data$POS))

      ## get context matrix
      mut_mat <-
        MutationalPatterns::mut_matrix(
          vcf_list = vcfs,
          ref_genome = pcgr_data[["assembly"]][["ref_genome"]],
          extension = 1)
      mut_mat <- mut_mat + 0.0001

      ## get reference signatures (COSMIC v3.2)
      all_reference_signatures <-
        MutationalPatterns::get_known_signatures(
        muttype = "snv",
        genome = stringr::str_replace(
          pcgr_data[["assembly"]][["grch_name"]], "grc", "GRC"
        ),
        incl_poss_artifacts =
          pcgr_config[["msigs"]][["include_artefact_signatures"]]
      )

      ## select subset of signatures based on those prevalent in tumor type/tissue
      selected_sigs <- intersect(
        colnames(all_reference_signatures),
        unique(prevalent_site_signatures$aetiology$signature_id)
      )
      selected_reference_signatures <-
        all_reference_signatures[, selected_sigs]

      ## reconstruct mutation profile from reference mutational signatures
      fit_ref <-
        MutationalPatterns::fit_to_signatures(mut_mat, selected_reference_signatures)

      ## assess how well the mutational profile can be reconstructed with the
      ## reference mutational signatures (cosine similarity)
      sim_original_reconstructed <-
        as.data.frame(
          MutationalPatterns::cos_sim_matrix(
            mut_mat, fit_ref[["reconstructed"]])) %>%
        magrittr::set_colnames("cosine_sim") %>%
        magrittr::set_rownames(NULL)

      ## assess the relative contribution of each reference mutational signature
      tot <- as.data.frame(
        stats::setNames(reshape2::melt(colSums(fit_ref[["contribution"]])),
                 c("tot"))) %>%
        dplyr::mutate(sample_id = as.character(rownames(.))) %>%
        magrittr::set_rownames(NULL)

      ## add information on aetiologies, and aggregate contributions
      ## pr. aetiology
      contributions_per_signature <-
        as.data.frame(stats::setNames(reshape2::melt(fit_ref[["contribution"]]),
                               c("signature_id", "sample_id",
                                 "contribution_raw"))) %>%
        dplyr::mutate(signature_id = as.character(.data$signature_id)) %>%
        dplyr::mutate(sample_id = as.character(.data$sample_id)) %>%
        dplyr::left_join(tot, by = "sample_id") %>%
        dplyr::mutate(prop_signature = round(as.numeric(.data$contribution_raw) / tot,
                                             digits = 3)) %>%
        dplyr::select(.data$signature_id, .data$sample_id, .data$prop_signature) %>%
        dplyr::filter(.data$prop_signature > 0) %>%
        dplyr::arrange(dplyr::desc(.data$prop_signature)) %>%
        dplyr::left_join(
          dplyr::select(pcgr_data[["mutational_signatures"]][["aetiologies"]],
                        .data$signature_id,
                        .data$aetiology, .data$comments, .data$aetiology_keyword),
          by = c("signature_id")) %>%
        dplyr::rename(group = .data$aetiology_keyword) %>%
        dplyr::mutate(
          contribution =
            paste0(round(.data$prop_signature * 100, digits = 2), "%")) %>%
        dplyr::distinct()

      contributions_per_group <- as.data.frame(
        contributions_per_signature %>%
          dplyr::group_by(.data$group) %>%
          dplyr::summarise(prop_group = sum(.data$prop_signature),
                         signature_id_group = paste(.data$signature_id, collapse=", "),
                         .groups = "drop")

      )
      ## FIX: if more than 18 aetiology types (unlikely), this will fail
      cols <- contributions_per_group %>%
        dplyr::arrange(dplyr::desc(.data$prop_group)) %>%
        dplyr::select(.data$group) %>%
        dplyr::distinct()

      color_vec <- utils::head(pcgrr::color_palette[["tier"]][["values"]], nrow(cols))

      names(color_vec) <- cols$group
      color_vec2 <- color_vec
      names(color_vec2) <- NULL
      cols <- cols %>% dplyr::mutate(col = color_vec2)
      contributions_per_signature <- contributions_per_signature %>%
        dplyr::left_join(cols, by = "group")

      contributions <- list()
      contributions[["per_group"]] <-  contributions_per_group
      contributions[["per_signature"]] <-  contributions_per_signature

      ## Get output for tab-separated file
      ## - contribution per signature id and reference signatures used
      if(!is.null(prevalent_site_signatures$aetiology) &
         NROW(contributions[["per_signature"]]) > 0){
        if("signature_id" %in% colnames(prevalent_site_signatures$aetiology)){
          reference_sigs <- paste(sort(prevalent_site_signatures$aetiology$signature_id),
                                  collapse=",")
          tsv_data <- contributions[["per_signature"]] %>%
            pcgrr::remove_cols_from_df(
              cnames = c("contribution","col","aetiology","comments")) %>%
            dplyr::mutate(
              all_reference_signatures = !type_specific,
              tumor_type = pcgr_config[["t_props"]][["tumor_type"]],
              reference_collection = "COSMIC_v32",
              reference_signatures = reference_sigs,
              fitting_accuracy =
                round(sim_original_reconstructed$cosine_sim * 100, digits = 1))
        }
      }

      vr <- vcfs[[sample_name]]
      GenomeInfoDb::seqlengths(vr) <-
        GenomeInfoDb::seqlengths(pcgr_data[["assembly"]][["bsg"]])[GenomeInfoDb::seqlevels(pcgr_data[["assembly"]][["bsg"]]) %in% unique(GenomeInfoDb::seqlevels(vr))]
      chromosomes <- utils::head(GenomeInfoDb::seqnames(pcgr_data[["assembly"]][["bsg"]]), 24)

      pcg_report_signatures[["result"]][["vr"]] <- vr
      pcg_report_signatures[["result"]][["mut_mat"]] <- mut_mat
      pcg_report_signatures[["result"]][["chromosomes"]] <- chromosomes
      pcg_report_signatures[["result"]][["contributions"]] <- contributions
      pcg_report_signatures[["result"]][["tsv"]] <- tsv_data
      pcg_report_signatures[["result"]][["reference_data"]] <- prevalent_site_signatures$aetiology
      pcg_report_signatures[["result"]][["scale_fill_values"]] <- color_vec
      pcg_report_signatures[["result"]][["scale_fill_names"]] <-
        names(color_vec)
      pcg_report_signatures[["result"]][["goodness_of_fit"]] <-
        round(sim_original_reconstructed$cosine_sim * 100, digits = 1)
    }else{
      pcg_report_signatures[["missing_data"]] <- TRUE
      if (length(vcfs[[1]]) > 0) {
        pcg_report_signatures[["variant_set"]][["all"]] <-
          as.data.frame(vcfs[[1]]) %>%
          dplyr::rename(POS = .data$start, CHROM = .data$seqnames) %>%
          dplyr::select(.data$CHROM, .data$POS, .data$REF, .data$ALT) %>%
          magrittr::set_rownames(NULL)
        log4r_info(
          paste0("Too few SNVs (n = ",
                 nrow(pcg_report_signatures[["variant_set"]][["all"]]),
                 ") for reconstruction of mutational signatures by ",
                 "MutationalPatterns, limit set to ",
                 pcgr_config[["msigs"]][["mutation_limit"]]))
      }
    }
  }

  return(pcg_report_signatures)
}


#' @export
get_prevalent_site_signatures <-
  function(site = "Any",
           custom_collection = NULL,
           pcgr_data = NULL,
           prevalence_pct = 5,
           incl_poss_artifacts = T) {

    if(is.null(custom_collection)){
      log4r_info(paste0(
        "Retrieving prevalent (prevalence >= ",
        prevalence_pct, " percent) reference signatures for ",
        site, ", using COSMIC v3.2 collection"))
    }
    log4r_info(paste0(
      "Inclusion of mutational signature artefacts (e.g. sequencing artefacts): ",
      incl_poss_artifacts))

    invisible(
      assertthat::assert_that(
        !is.null(pcgr_data[["mutational_signatures"]][["aetiologies"]]),
        msg =
          "Cannot load ref. aetiologies (COSMIC v3.2) of mutational signatures"))
    invisible(
      assertthat::assert_that(
        is.data.frame(pcgr_data[["mutational_signatures"]][["aetiologies"]]),
        msg = "Reference aetiologies must be of type data.frame()"))
    invisible(
      assertthat::assert_that(
        prevalence_pct == 0 |
          prevalence_pct == 2 | prevalence_pct == 5 |
          prevalence_pct == 10 | prevalence_pct == 15 |
          prevalence_pct == 20,
        msg = "Argument 'prevalence_pct' must be any of '0, 2, 5, 10, 15 or 20'"))

    valid_signature_ids <-
      unique(pcgr_data[["mutational_signatures"]][["aetiologies"]]$signature_id)
    signatures_prevalence <- data.frame()

    if(!is.null(custom_collection)){
      invisible(
        assertthat::assert_that(
          is.character(custom_collection),
          msg = "Argument 'custom_collection' must be a character vector"))

      log4r_info(paste0(
        "Retrieving reference signatures from COSMIC v3.2 collection based on user-defined collection (",
        paste(unique(custom_collection), collapse=", "), ")")
      )
      i <- 1
      while(i <= length(custom_collection)){
        if(!(custom_collection[i] %in% valid_signature_ids)){
          log4r_warn(paste0("Could not find specified custom signature id  '",
                                    custom_collection[i], "' in COSMIC v3.2 reference collection",
                                    " - ignoring"))
        }
        i <- i + 1
      }

      signatures_prevalence <-
        pcgr_data[["mutational_signatures"]][["aetiologies"]] %>%
        dplyr::select(.data$signature_id,
                      .data$aetiology_keyword,
                      .data$aetiology,
                      .data$associated_signatures,
                      .data$comments) %>%
        dplyr::filter(.data$signature_id %in% custom_collection) %>%
        dplyr::distinct()

    }else{

      unique_sites_with_signature_prevalence <-
        unique(pcgr_data[["mutational_signatures"]][["aetiologies"]][["primary_site"]])
      if (!(site %in% unique_sites_with_signature_prevalence)) {
        log4r_info(
          paste0("Primary tumor site '", site, "' ",
                 "does not have any signatures with significant ",
                 "prevalence - considering all"))
        signatures_prevalence <-
          pcgr_data[["mutational_signatures"]][["aetiologies"]] %>%
          dplyr::select(.data$signature_id, .data$aetiology_keyword, .data$aetiology,
                        .data$associated_signatures, .data$comments) %>%
          dplyr::distinct()
      }else{
        signatures_prevalence <-
          pcgr_data[["mutational_signatures"]][["aetiologies"]] %>%
          dplyr::filter(.data$primary_site == site) %>%
          dplyr::select(.data$signature_id,
                        .data$primary_site,
                        .data$prevalence_pct,
                        .data$prevalence_above_5pct,
                        .data$prevalence_above_10pct,
                        .data$prevalence_above_15pct,
                        .data$prevalence_above_20pct,
                        .data$aetiology_keyword, .data$aetiology,
                        .data$associated_signatures, .data$comments) %>%
          dplyr::distinct()

        if (prevalence_pct > 0) {
          if (prevalence_pct == 5) {
            signatures_prevalence <- signatures_prevalence %>%
              dplyr::filter(.data$prevalence_above_5pct == T |
                              is.na(.data$prevalence_above_5pct))
          }else if (prevalence_pct == 10) {
            signatures_prevalence <- signatures_prevalence %>%
              dplyr::filter(.data$prevalence_above_10pct == T |
                              is.na(.data$prevalence_above_10pct))
          }
          else if (prevalence_pct == 15) {
            signatures_prevalence <- signatures_prevalence %>%
              dplyr::filter(.data$prevalence_above_15pct == T |
                              is.na(.data$prevalence_above_15pct))
          }else if (prevalence_pct == 20) {
            signatures_prevalence <- signatures_prevalence %>%
              dplyr::filter(.data$prevalence_above_20pct == T |
                              is.na(.data$prevalence_above_20pct))
          }
        }
        signatures_prevalence <- signatures_prevalence %>%
          dplyr::select(-c(.data$primary_site, .data$prevalence_above_5pct,
                           .data$prevalence_above_10pct, .data$prevalence_above_15pct,
                           .data$prevalence_above_20pct)) %>%
          dplyr::distinct() %>%
          dplyr::arrange(dplyr::desc(.data$prevalence_pct)) %>%
          dplyr::select(-.data$prevalence_pct)
      }
    }

    if(incl_poss_artifacts == F){
      signatures_prevalence <- signatures_prevalence %>%
        dplyr::filter(!stringr::str_detect(.data$aetiology_keyword,"artefact"))
    }
    signatures_prevalence <- signatures_prevalence %>%
      dplyr::distinct()

    ## Subset signature matrix - keeping only columns (signatures)
    ## to those defined by primary site/custom collection
    sigs <- unique(signatures_prevalence$signature_id)
    log4r_info(paste0("Limiting reference collection to signatures: ",
                              paste(sigs, collapse = ", ")))

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
generate_report_data_rainfall <- function(variant_set, colors = NULL,
                                          autosomes = F, build = "grch37") {

  pcg_report_rainfall <- pcgrr::init_report(class = "rainfall")
  if(NROW(variant_set) == 0){
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

  log4r_info("------")
  log4r_info(paste0("Calculating data for rainfall plot"))


  sbs_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(6, "Dark2")
    # colors <- c(
    #   "#2EBAED", "#000000", "#DE1C14",
    #   "#D4D2D2", "#ADCC54", "#F0D0CE")
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
    variant_set <- variant_set %>%
      dplyr::mutate(CHROM = paste0("chr", .data$CHROM))
  }

  dat <- dplyr::select(variant_set, .data$CHROM, .data$POS,
                       .data$REF, .data$ALT, .data$VARIANT_CLASS) %>%
    dplyr::filter(.data$VARIANT_CLASS == "SNV")

  if (nrow(dat) < 10 | length(chromosome_names) < 2) {
    log4r_info(
      paste0("Too few variants (< 10) and chromosomes ",
             " represented (< 2) - skipping rainfall plot"))
    pcg_report_rainfall[["eval"]] <- F
  }else{
    dat <- dat %>%
      pcgrr::assign_mutation_type() %>%
      dplyr::mutate(
        MUTATION_TYPE =
          dplyr::if_else(stringr::str_detect(.data$MUTATION_TYPE, "^C>"),
                         stringr::str_replace(.data$MUTATION_TYPE,
                                              ":[A-Z]>[A-Z]$", ""),
                         as.character(.data$MUTATION_TYPE))) %>%
      dplyr::mutate(
        MUTATION_TYPE =
          dplyr::if_else(stringr::str_detect(.data$MUTATION_TYPE, "^A>"),
                         stringr::str_replace(.data$MUTATION_TYPE,
                                              "^[A-Z]>[A-Z]:", ""),
                         as.character(.data$MUTATION_TYPE))) %>%
      pcgrr::sort_chromosomal_segments()

    bsg <- BSgenome.Hsapiens.UCSC.hg19
    chr_length <- utils::head(GenomeInfoDb::seqlengths(bsg), 24)
    chromosomes <- utils::head(GenomeInfoDb::seqnames(bsg), 24)
    if (build == "grch38") {
      bsg <- BSgenome.Hsapiens.UCSC.hg38
      chr_length <- utils::head(GenomeInfoDb::seqlengths(bsg), 24)
    }
    if (autosomes == T) {
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
