#' Function that retrieves relative estimates of known somatic signatures
#' from a single tumor
#'
#' @param mut_data data frame with somatic mutations (VCF_SAMPLE_ID, CHROM,
#' POS, REF, ALT)
#' @param sample_name sample name
#' @param normalization_method metod for normalization of context
#' counts (deconstructSigs)
#' @param cosmic_cancertypes_aetiologies list of known signaturea and
#' associated etiologies/cancertypes
#' @param signature_limit max number of contributing signatures
#' @param associated_signatures limit search spaced vector of signatures
#' that are  c("Signature.29","Signature.23")
#' @param signature_cutoff discard any signature contributions with weight
#' less than this amount
#' @param bsg genome sequence object (BSgenome.Hsapiens.UCSC)
#' @param weight_precision number of significant digits in signature weight
#'
#'
signature_contributions_single_sample <-
  function(sample_calls,
           sample_name,
           normalization_method = "default",
           cosmic_signatures_aetiologies = NULL,
           signature_limit = 6,
           associated_signatures = NULL,
           signature_cutoff = 0.06,
           bsgenome = NULL,
           weight_precision = 3) {

  assertable::assert_colnames(
    sample_calls, c("VCF_SAMPLE_ID", "CHROM", "POS", "REF", "ALT"),
    only_colnames = F, quiet = T)
  sample_calls <- pcgrr::get_valid_chromosomes(sample_calls)
  sigs.input <- deconstructSigs::mut.to.sigs.input(
    mut.ref = sample_calls,
    sample.id = "VCF_SAMPLE_ID",
    chr = "CHROM",
    pos = "POS",
    ref = "REF",
    alt = "ALT", bsg = bsgenome)
  all_signatures <- paste0("Signature.", rep(1:30))
  if (!is.null(associated_signatures)) {
    all_signatures <- all_signatures
  }else{
    all_signatures <- associated_signatures
  }

  sample_1 <-
    deconstructSigs::whichSignatures(tumor.ref = sigs.input,
                                     sample.id = sample_name,
                                     associated = all_signatures,
                                     signature.cutoff = signature_cutoff,
                                     signatures.limit = signature_limit,
                                     signatures.ref = signatures.cosmic,
                                     contexts.needed = T,
                                     tri.counts.method = normalization_method)

  nonzero_signatures <-
    sample_1$weights[which(colSums(sample_1$weights != 0) > 0)]
  n <- 1
  signature_contributions <- NULL
  while (n <= ncol(nonzero_signatures)) {
    df <- data.frame("sample_name" = sample_name,
                     "signature_id" =
                       stringr::str_replace(colnames(nonzero_signatures)[n],
                                            "ignature\\.", ""),
                     "weight" = as.numeric(nonzero_signatures[, n]))
    signature_contributions <- rbind(signature_contributions, df)
    rlogging::message(paste0("Inferred weighted contribution of ",
                             df$signature_id, ": ",
                             round(df$weight,
                                   digits = weight_precision)))
    n <- n + 1
  }

  signature_contributions <-
    rbind(signature_contributions, data.frame("sample_name" = sample_name,
                                              "signature_id" = "unknown",
                                              "weight" = sample_1$unknown))
  signature_rows <- as.numeric(
    stringr::str_replace(
      as.character(
        signature_contributions[signature_contributions$signature_id != "unknown", ]$signature_id), "S", ""))
  weight_df <-
    data.frame("Signature_ID" =
                 as.character(signature_contributions$signature_id),
               "Weight" =
                 round(as.numeric(signature_contributions$weight),
                                digits = 3),
               stringsAsFactors = F)

  signatures_cancertypes_aetiologies <-
    cosmic_signatures_aetiologies[signature_rows, ] %>%
    dplyr::left_join(weight_df, by = c("Signature_ID")) %>%
    dplyr::arrange(desc(Weight)) %>%
    dplyr::select(Signature_ID, Weight, Cancer_types,
                  Proposed_aetiology, Comments, Keyword) %>%
    dplyr::mutate(Trimer_normalization_method = normalization_method,
                  Sample_Name = sample_name)

  if (nrow(
    signature_contributions[signature_contributions$signature_id == "unknown", ]) > 0) {
    unknown_df <- as.data.frame(
      signature_contributions %>%
        dplyr::filter(signature_id == "unknown") %>%
        dplyr::rename(Weight = weight) %>%
        dplyr::mutate(Weight = round(as.numeric(Weight),
                                     digits = weight_precision)) %>%
        dplyr::mutate(Sample_Name = as.character(sample_name),
                      Trimer_normalization_method =
                        as.character(normalization_method)) %>%
        dplyr::mutate(Signature_ID = as.character(signature_id)) %>%
        dplyr::select(-c(signature_id, sample_name))
    )
    signatures_cancertypes_aetiologies <-
      dplyr::bind_rows(signatures_cancertypes_aetiologies, unknown_df)
  }
  signatures_cancertypes_aetiologies <-
    signatures_cancertypes_aetiologies %>%
    dplyr::arrange(desc(Weight))

  return(list(
    deconstructsigs_which_signatures = sample_1,
    cancertypes_aetiologies = signatures_cancertypes_aetiologies))
}

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
generate_report_data_signatures_mp <-
  function(vcf_fname, pcgr_data, sample_name,
           pcgr_config,
           type_specific = T) {

  rlogging::message("------")
  rlogging::message("Identifying weighted contributions of reference ",
                    "mutational signatures (COSMIC v3) using ",
                    "MutationalPatterns")
  assay <- tolower(pcgr_config$assay_props$type)
  pcg_report_signatures <-
    pcgrr::init_report(pcgr_config, sample_name, class = "m_signature_mp")
  reference_data <-
    pcgrr::get_reference_signature_data(
      pcgr_config[["t_props"]][["tumor_type"]],
      assay = assay)
  if(type_specific == F){
    reference_data <-
      pcgrr::get_reference_signature_data("Any", assay = assay)
  }
  if(file.exists(vcf_fname)){
    vcfs <- suppressWarnings(
      MutationalPatterns::read_vcfs_as_granges(
        vcf_fname, sample_name,
        pcgr_data[["assembly"]][["ref_genome"]]))

    rlogging::message(paste0("Number of SNVs for signature analysis: ",
                             length(vcfs[[1]])))

    pcg_report_signatures[["eval"]] <- TRUE

    if (length(vcfs[[1]]) >= pcgr_config[["msigs"]][["mutation_limit"]]) {
      pcg_report_signatures[["variant_set"]][["all"]] <-
        as.data.frame(vcfs[[1]]) %>%
        dplyr::rename(POS = start, CHROM = seqnames) %>%
        dplyr::select(CHROM, POS, REF, ALT) %>%
        magrittr::set_rownames(NULL)

      ## get context matrix
      mut_mat <-
        MutationalPatterns::mut_matrix(
          vcf_list = vcfs,
          ref_genome = pcgr_data[["assembly"]][["ref_genome"]])
      mut_mat <- mut_mat + 0.0001

      ## reconstruct mutation profile from reference mutational signatures
      fit_ref <- MutationalPatterns::fit_to_signatures(
        mut_mat, reference_data$matrix)

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
        setNames(reshape2::melt(colSums(fit_ref[["contribution"]])),
                 c("tot"))) %>%
        dplyr::mutate(sample_id = as.character(rownames(.))) %>%
        magrittr::set_rownames(NULL)

      ## add information on aetiologies, and aggregate contributions
      ## pr. aetiology
      contributions_per_signature <-
        as.data.frame(setNames(reshape2::melt(fit_ref[["contribution"]]),
                               c("signature_id", "sample_id",
                                 "contribution_raw"))) %>%
        dplyr::mutate(signature_id = as.character(signature_id)) %>%
        dplyr::mutate(sample_id = as.character(sample_id)) %>%
        dplyr::left_join(tot, by = "sample_id") %>%
        dplyr::mutate(prop_signature = round(as.numeric(contribution_raw) / tot,
                                             digits = 3)) %>%
        dplyr::select(signature_id, sample_id, prop_signature) %>%
        dplyr::filter(prop_signature > 0) %>%
        dplyr::arrange(desc(prop_signature)) %>%
        dplyr::left_join(
          dplyr::select(pcgr_data[["mutational_signatures"]][["aetiologies_67"]],
                        signature_id,
                        aetiology, comments, aetiology_keyword),
          by = c("signature_id")) %>%
        dplyr::rename(group = aetiology_keyword) %>%
        dplyr::mutate(
          contribution =
            paste0(round(prop_signature * 100, digits = 2), "%")) %>%
        dplyr::distinct()

      contributions_per_group <- as.data.frame(
        contributions_per_signature %>%
          dplyr::group_by(group) %>%
          dplyr::summarise(prop_group = sum(prop_signature),
                         signature_id_group = paste(signature_id, collapse=", "),
                         .groups = "drop")

      )
      ## FIX: if more than 18 aetiology types (unlikely), this will fail
      cols <- contributions_per_group %>%
        dplyr::arrange(desc(prop_group)) %>%
        dplyr::select(group) %>%
        dplyr::distinct()

      color_vec <- head(pcgrr::color_palette[["tier"]][["values"]], nrow(cols))

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
      if(!is.null(reference_data$aetiology) &
         NROW(contributions[["per_signature"]]) > 0){
        if("signature_id" %in% colnames(reference_data$aetiology)){
          reference_sigs <- paste(sort(reference_data$aetiology$signature_id),
                                  collapse=",")
          tsv_data <- contributions[["per_signature"]] %>%
            pcgrr::remove_cols_from_df(
              cnames = c("contribution","col","aetiology","comments")) %>%
            dplyr::mutate(
              all_reference_signatures = !type_specific,
              tumor_type = pcgr_config[["t_props"]][["tumor_type"]],
              reference_collection = paste0("COSMIC_v3_",toupper(assay)),
              reference_signatures = reference_sigs,
              fitting_accuracy =
                round(sim_original_reconstructed$cosine_sim * 100, digits = 1))
        }
      }

      vr <- vcfs[[sample_name]]
      GenomeInfoDb::seqlengths(vr) <-
        GenomeInfoDb::seqlengths(pcgr_data[["assembly"]][["bsg"]])[GenomeInfoDb::seqlevels(pcgr_data[["assembly"]][["bsg"]]) %in% unique(GenomeInfoDb::seqlevels(vr))]
      chromosomes <- head(GenomeInfoDb::seqnames(pcgr_data[["assembly"]][["bsg"]]), 24)

      pcg_report_signatures[["result"]][["vr"]] <- vr
      pcg_report_signatures[["result"]][["mut_mat"]] <- mut_mat
      pcg_report_signatures[["result"]][["chromosomes"]] <- chromosomes
      pcg_report_signatures[["result"]][["contributions"]] <- contributions
      pcg_report_signatures[["result"]][["tsv"]] <- tsv_data
      pcg_report_signatures[["result"]][["reference_data"]] <- reference_data
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
          dplyr::rename(POS = start, CHROM = seqnames) %>%
          dplyr::select(CHROM, POS, REF, ALT) %>%
          magrittr::set_rownames(NULL)
        rlogging::message(
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


get_reference_signature_data <- function(site, assay = "wes",
                                         prevalence_pct = 2,
                                         ignore_seq_artefacts = T) {

  rlogging::message("Retrieving prevalent (prevalence >= ",
                    prevalence_pct, " percent) reference signatures for ",
                    site, ", using COSMIC v3 collection (", toupper(assay), ")")
  #assertthat::assert_that(!)
  invisible(
    assertthat::assert_that(
      !is.null(pcgr_data[["mutational_signatures"]][["aetiologies_67"]]),
      msg =
        "Cannot load ref. collection (COSMIC v3) of mutational signatures"))
  invisible(
    assertthat::assert_that(
      is.data.frame(pcgr_data[["mutational_signatures"]][["aetiologies_67"]]),
      msg = "Reference collection must be of type data.frame()"))
  invisible(
    assertthat::assert_that(
      assay == "wes" | assay == "wgs",
      msg = "Argument 'assay' must be either 'wes' or 'wgs'"))
  invisible(
    assertthat::assert_that(
      prevalence_pct == 2 | prevalence_pct == 5 |
        prevalence_pct == 10 | prevalence_pct == 15 |
        prevalence_pct == 20,
      msg = "Argument 'prevalence_pct' must be any of '2, 5, 10, 15 or 20'"))

  ##probabilities of mutation contexts in each of 67 cancer signatures
  ## (from COSMIC/ICGC)
  cancer_sigs <- pcgr_data[["mutational_signatures"]][["probabilities_wes"]]
  if (assay == "wgs") {
    cancer_sigs <- pcgr_data[["mutational_signatures"]][["probabilities_wgs"]]
  }
  rownames(cancer_sigs) <- paste0(substr(cancer_sigs$SubType, 1, 1),
                                  "[", cancer_sigs$Type, "]",
                                  substr(cancer_sigs$SubType, 3, 3))
  cancer_sigs$SubType <- NULL
  cancer_sigs$Type <- NULL
  signatures_prevalence <- data.frame()

  unique_sites_with_signature_prevalence <-
    unique(pcgr_data[["mutational_signatures"]][["aetiologies_67"]][["primary_site"]])
  if (!(site %in% unique_sites_with_signature_prevalence)) {
    rlogging::message(
      paste0("Primary tumor site '", site, "'",
             "does not have any signatures with significant ",
             "prevalence - considering all"))
    signatures_prevalence <-
      pcgr_data[["mutational_signatures"]][["aetiologies_67"]] %>%
      dplyr::select(signature_id, aetiology_keyword, aetiology,
                    associated_signatures, comments) %>%
      dplyr::filter(signature_id != "SBS84" & signature_id != "SBS85") %>%
      dplyr::distinct()
  }else{
    signatures_prevalence <-
      pcgr_data[["mutational_signatures"]][["aetiologies_67"]] %>%
      dplyr::filter(primary_site == site) %>%
      dplyr::select(signature_id, primary_site, prevalence_pct,
                    prevalence_above_5pct,
                    prevalence_above_10pct, prevalence_above_15pct,
                    prevalence_above_20pct, aetiology_keyword, aetiology,
                    associated_signatures, comments) %>%
      dplyr::distinct()

    if (prevalence_pct > 0) {
      if (prevalence_pct == 5) {
        signatures_prevalence <- signatures_prevalence %>%
          dplyr::filter(prevalence_above_5pct == T)
      }else if (prevalence_pct == 10) {
        signatures_prevalence <- signatures_prevalence %>%
          dplyr::filter(prevalence_above_10pct == T)
      }
      else if (prevalence_pct == 15) {
        signatures_prevalence <- signatures_prevalence %>%
          dplyr::filter(prevalence_above_15pct == T)
      }else if (prevalence_pct == 20) {
        signatures_prevalence <- signatures_prevalence %>%
          dplyr::filter(prevalence_above_20pct == T)
      }
    }
    signatures_prevalence <- signatures_prevalence %>%
      dplyr::select(-c(primary_site, prevalence_above_5pct,
                       prevalence_above_10pct, prevalence_above_15pct,
                       prevalence_above_20pct)) %>%
      dplyr::distinct() %>%
      dplyr::filter(signature_id != "SBS84" &
                      signature_id != "SBS85") %>%
      dplyr::arrange(desc(prevalence_pct)) %>%
      dplyr::select(-prevalence_pct)
  }

  sigs <- unique(signatures_prevalence$signature_id)
  signatures_prevalence <- signatures_prevalence %>%
    dplyr::distinct()
  rlogging::message(paste0("Limiting reference collection to signatures: ",
                           paste(sigs, collapse = ", ")))
  cancer_sigs <- dplyr::select(cancer_sigs, dplyr::one_of(sigs))


  result <- list("matrix" = as.matrix(cancer_sigs),
                 "aetiology" = signatures_prevalence)

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
generate_report_data_rainfall <- function(variant_set, colors = NULL,
                                          autosomes = F, build = "grch37") {

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

  rlogging::message("------")
  rlogging::message(paste0("Calculating data for rainfall plot"))

  pcg_report_rainfall <- pcgrr::init_report(class = "rainfall")


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
      dplyr::mutate(CHROM = paste0("chr", CHROM))
  }

  dat <- dplyr::select(variant_set, CHROM, POS,
                       REF, ALT, VARIANT_CLASS) %>%
    dplyr::filter(VARIANT_CLASS == "SNV")

  if (nrow(dat) < 10 | length(chromosome_names) < 2) {
    rlogging::message("Too few variants (< 10) and chromosomes ",
                      "represented (< 2) - skipping rainfall plot")
    pcg_report_rainfall[["eval"]] <- F
  }else{
    dat <- dat %>%
      pcgrr::assign_mutation_type() %>%
      dplyr::mutate(
        MUTATION_TYPE =
          dplyr::if_else(stringr::str_detect(MUTATION_TYPE, "^C>"),
                         stringr::str_replace(MUTATION_TYPE,
                                              ":[A-Z]>[A-Z]$", ""),
                         as.character(MUTATION_TYPE))) %>%
      dplyr::mutate(
        MUTATION_TYPE =
          dplyr::if_else(stringr::str_detect(MUTATION_TYPE, "^A>"),
                         stringr::str_replace(MUTATION_TYPE,
                                              "^[A-Z]>[A-Z]:", ""),
                         as.character(MUTATION_TYPE))) %>%
      pcgrr::sort_chromosomal_segments()

    bsg <- BSgenome.Hsapiens.UCSC.hg19
    chr_length <- head(GenomeInfoDb::seqlengths(bsg), 24)
    chromosomes <- head(GenomeInfoDb::seqnames(bsg), 24)
    if (build == "grch38") {
      bsg <- BSgenome.Hsapiens.UCSC.hg38
      chr_length <- head(GenomeInfoDb::seqlengths(bsg), 24)
    }
    if (autosomes == T) {
      chr_length <- head(chr_length, 22)
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
      chr_subset <- dplyr::filter(dat, CHROM == chromosomes[i])
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
